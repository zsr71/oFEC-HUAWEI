#include "newcode/ofec_decoder.hpp"
#include "newcode/chase256.hpp" // 保留 Chase 头；本文档内有三参前向声明
#include "newcode/bch_255_239.hpp"
#include "newcode/ofec_decoder_hard.hpp"
#include "newcode/decoder_core.hpp"
#include "newcode/llr_utils.hpp"
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <type_traits>
#include <limits>

namespace newcode {

// 三参版本前向声明（定义在 chase256.cpp，已对常用类型显式实例化）
template<typename LLR>
void chase_decode_256(const LLR* Lin256, const LLR* Lch256, LLR* Y2_256, const Params& p);

template <typename LLR>
static bool tile_should_early_stop(const Matrix<LLR>& lin_matrix)
{
  // 预期列数为 256 = 128(旧) + 111(新) + 16(BCH校验) + 1(整体奇偶)
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  if (cols < 256) return false; // 保守：尺寸异常则不早停

  std::array<uint8_t, 255> hard255;
  std::array<uint8_t, 255> decoded255;

  for (size_t r = 0; r < rows; ++r)
  {
    // 1) LLR -> 硬判到 {0,1}，其中 1 表示比特“1”
    for (int j = 0; j < 255; ++j) {
      const float v = llr_to_float(lin_matrix[r][static_cast<size_t>(j)]);
      hard255[static_cast<size_t>(j)] = (v < 0.0f) ? 1u : 0u;
    }

    // 2) 判定是否为合法 BCH(255,239) 码字
    if (!bch_255_239_decode_hiho_cw_255(hard255.data(), decoded255.data()))
      return false; // 该行无法成功译码，不能早停

    // 3) 扩展整体偶校验（第 256 位在索引 255）
    uint8_t parity255 = 0u;
    for (int j = 0; j < 255; ++j) parity255 ^= (hard255[static_cast<size_t>(j)] & 1u);
    const uint8_t overall = (llr_to_float(lin_matrix[r][255]) < 0.0f) ? 1u : 0u;

    // 扩展码要求 255 位与整体位异或为 0（偶校验通过）
    if ((parity255 ^ overall) != 0u)
      return false; // 整体校验失败，不能早停
  }

  // 所有行均满足：BCH 综合全零且整体偶校验通过 => 触发早停
  return true;
}


template <typename LLR>
TileProcessResult<LLR> process_tile(const Matrix<LLR>& tile_in,
                                    const Matrix<LLR>& ch_tile,
                                    const Params& p,
                                    size_t tile_top_row_global,
                                    bool use_hard_decode)
{
  constexpr int B         = static_cast<int>(Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(Params::NUM_SUBBLOCK_COLS * B);            // 128
  constexpr int K         = static_cast<int>(Params::BCH_K);                             // 239
  constexpr int TAKE_BITS = K - N;                                                       // 111
  constexpr int BCH_PAR   = static_cast<int>(Params::BCH_PARITY_BITS);                   // 16
  constexpr int OVR_IDX   = static_cast<int>(Params::BCH_OVERALL_IDX);                   // 255

  const size_t H = tile_in.rows();
  const size_t W = tile_in.cols();
  assert(W == static_cast<size_t>(N));
  assert(ch_tile.rows() == H && ch_tile.cols() == W);

  Matrix<LLR> tile_out = tile_in; // 回写外信息时覆盖相关位置

  const int SBR = p.CHASE_SBR;
  if (SBR != 1 && SBR != 2)
      throw std::invalid_argument("process_tile: CHASE_SBR must be 1 or 2.");

  const size_t rows_to_decode = static_cast<size_t>(SBR) * static_cast<size_t>(B);
  if (rows_to_decode == 0) {
      return TileProcessResult<LLR>{tile_out, false};
  }
  const size_t decoder_cols = static_cast<size_t>(2 * N);
  Matrix<LLR> lin_matrix(rows_to_decode, decoder_cols);
  Matrix<LLR> lch_matrix(rows_to_decode, decoder_cols);
  std::vector<size_t> row_local_lookup(rows_to_decode, 0);
  std::vector<size_t> row_global_lookup(rows_to_decode, 0);

  // 遍历底部 SBR 个 sub-block rows：从最底开始，组装 decoder 输入矩阵
  for (int s = 0; s < SBR; ++s)
  {
      const size_t sbr_row0_local = H - static_cast<size_t>((s + 1) * B);
      for (int r_off = 0; r_off < B; ++r_off)
      {
          const size_t row_idx   = static_cast<size_t>(s * B + r_off);
          const size_t row_local = sbr_row0_local + static_cast<size_t>(r_off);
          const size_t row_global = tile_top_row_global + row_local;

          row_local_lookup[row_idx]  = row_local;
          row_global_lookup[row_idx] = row_global;

          const long R = static_cast<long>(row_global / static_cast<size_t>(B));
          const int  r = static_cast<int>(row_global % static_cast<size_t>(B));

          // (A) 0..N-1，"旧信息" 128 位（跨块映射）
          for (int k = 0; k < N; ++k)
          {
              const long br = (R ^ 1L) - static_cast<long>(2 * p.NUM_GUARD_SUBROWS)
                              - static_cast<long>(2 * (N / B))
                              + static_cast<long>(2 * (k / B));
              const long bc = static_cast<long>(k / B);
              const long bit_row_in_block = static_cast<long>((k % B) ^ r);
              const long bit_col_in_block = static_cast<long>(r);

              const long rr_global = br * B + bit_row_in_block;
              const long cc_global = bc * B + bit_col_in_block;

              const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
              const long cc_local2 = cc_global; // 0..N-1

              if (rr_local2 >= 0 && rr_local2 < static_cast<long>(H) &&
                  cc_local2 >= 0 && cc_local2 < static_cast<long>(W))
              {
                  const float Lch = llr_to_float(ch_tile[static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)]);
                  const float La  = llr_to_float(tile_in [static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)]);
                  lin_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch + La);
                  lch_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch);
              }
              else {
                  throw std::out_of_range("process_tile: old info position out of tile range.");
              }
          }

          // (B) N..K-1，111 个"新信息"（当前行）
          for (int i = 0; i < TAKE_BITS; ++i) {
              const float Lch = llr_to_float(ch_tile[row_local][static_cast<size_t>(i)]);
              const float La  = llr_to_float(tile_in [row_local][static_cast<size_t>(i)]);
              const size_t col = static_cast<size_t>(N + i);
              lin_matrix[row_idx][col] = llr_from_float<LLR>(Lch + La);
              lch_matrix[row_idx][col] = llr_from_float<LLR>(Lch);
          }

          // (C) 239..254，16 个 BCH 校验位（当前行）
          for (int j = 0; j < BCH_PAR; ++j) {
              const size_t col = static_cast<size_t>(K + j);
              const float Lch = llr_to_float(ch_tile[row_local][static_cast<size_t>(TAKE_BITS + j)]);
              const float La  = llr_to_float(tile_in [row_local][static_cast<size_t>(TAKE_BITS + j)]);
              lin_matrix[row_idx][col] = llr_from_float<LLR>(Lch + La);
              lch_matrix[row_idx][col] = llr_from_float<LLR>(Lch);
          }

          // (D) 255，整体偶校验位（当前行）
          {
              const size_t col = static_cast<size_t>(OVR_IDX);
              const size_t src_col = static_cast<size_t>(TAKE_BITS + BCH_PAR);
              const float Lch = llr_to_float(ch_tile[row_local][src_col]);
              const float La  = llr_to_float(tile_in [row_local][src_col]);
              lin_matrix[row_idx][col] = llr_from_float<LLR>(Lch + La);
              lch_matrix[row_idx][col] = llr_from_float<LLR>(Lch);
          }
      } // r_off
  } // s

  bool early_stop_triggered = tile_should_early_stop(lin_matrix);

  auto decoder_res = Decoder_Core(lin_matrix, lch_matrix, use_hard_decode, p);

  // ====== 回写外信息到 tile_out（作为下一轮/下一组件先验） ======
  for (int s = 0; s < SBR; ++s)
  {
      for (int r_off = 0; r_off < B; ++r_off)
      {
          const size_t row_idx = static_cast<size_t>(s * B + r_off);
          if (row_idx >= decoder_res.produced_rows.size())
              continue;
          if (!decoder_res.produced_rows[row_idx])
              continue; // 仅选择硬判时且失败，保留原外信息

          const size_t row_local  = row_local_lookup[row_idx];
          const size_t row_global = row_global_lookup[row_idx];
          const auto&  lout_row   = decoder_res.lout[row_idx];

          // 当前行：111 新信息 + 16 BCH + 1 overall
          for (int i = 0; i < TAKE_BITS; ++i)
              tile_out[row_local][static_cast<size_t>(i)] = lout_row[static_cast<size_t>(N + i)];
          for (int j = 0; j < BCH_PAR; ++j)
              tile_out[row_local][static_cast<size_t>(TAKE_BITS + j)] = lout_row[static_cast<size_t>(K + j)];
          tile_out[row_local][static_cast<size_t>(TAKE_BITS + BCH_PAR)] = lout_row[static_cast<size_t>(OVR_IDX)];

          // 旧信息 128：把 Y2[0..N-1] 按映射散布回对应 (rr,cc)
          const long R = static_cast<long>(row_global / static_cast<size_t>(B));
          const int  r = static_cast<int>(row_global % static_cast<size_t>(B));

          for (int k = 0; k < N; ++k)
          {
              const long br = (R ^ 1L) - static_cast<long>(2 * p.NUM_GUARD_SUBROWS)
                              - static_cast<long>(2 * (N / B))
                              + static_cast<long>(2 * (k / B));
              const long bc = static_cast<long>(k / B);
              const long bit_row_in_block = static_cast<long>((k % B) ^ r);
              const long bit_col_in_block = static_cast<long>(r);

              const long rr_global = br * B + bit_row_in_block;
              const long cc_global = bc * B + bit_col_in_block;

              const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
              const long cc_local2 = cc_global;

              const bool in_range =
                  (rr_local2 >= 0 && rr_local2 < static_cast<long>(H) &&
                  cc_local2 >= 0 && cc_local2 < static_cast<long>(W));

              // —— Debug：严格要求必须命中；Release：可选择抛异常或计数 —— //
              assert(in_range && "process_tile: write-back out of tile range");

              tile_out[static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)] =
                      lout_row[static_cast<size_t>(k)];
          }
      } // r_off
  } // s

  return TileProcessResult<LLR>{std::move(tile_out), early_stop_triggered};
} // process_tile

  // ========== 窗口处理 ==========
template <typename LLR>
void process_window(Matrix<LLR>& work_llr,
                    const Matrix<LLR>& channel_llr,
                    size_t win_start, size_t win_end, const Params& p,
                    size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                    std::vector<TileEarlyStopCounter>* tile_stats)
{
  const size_t N = Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;

  for (size_t t = 0; t < TILES_PER_WIN; ++t)
  {
        // 让 t=0 对应窗口最底部的 tile
        const size_t tile_bottom_row = win_end  - t * tile_height_rows;
        const size_t tile_top_row    = tile_bottom_row + 1 - tile_height_rows;


        const size_t tile_height_rows = tile_bottom_row - tile_top_row + 1;
        Matrix<LLR> tile_in(tile_height_rows, N);
        Matrix<LLR> ch_tile(tile_height_rows, N);

        for (size_t r = 0; r < tile_height_rows; ++r) {
            for (size_t c = 0; c < N; ++c) {
                tile_in[r][c]  = work_llr[tile_top_row + r][c];
                ch_tile[r][c]  = channel_llr[tile_top_row + r][c];
            }
        }

        // 参数选取
        auto pick_float = [](const std::vector<float>& tbl, size_t idx, float fallback) -> float {
            return (idx < tbl.size()) ? tbl[idx] : fallback;
        };
        auto pick_int = [](const std::vector<int>& tbl, size_t idx, int fallback) -> int {
            return (idx < tbl.size()) ? tbl[idx] : fallback;
        };

        Params tile_params = p;
        tile_params.beta = pick_float(p.beta_list, t, p.beta);
        tile_params.ALPHA = pick_float(p.ALPHA_LIST, t, p.ALPHA);

    const bool use_hard = pick_int(p.HARD_TILE_LIST, t, p.HARD_DECODE_DEFAULT ? 1 : 0) != 0;

    TileProcessResult<LLR> tile_result = process_tile<LLR>(tile_in, ch_tile, tile_params,
                                                           /*tile_top_row_global=*/tile_top_row,
                                                           /*use_hard_decode=*/use_hard);

    if (tile_stats && t < tile_stats->size()) {
      auto& counter = (*tile_stats)[t];
      counter.total += 1;
      if (tile_result.early_stop_triggered) {
        counter.triggered += 1;
      }
    }

    for (size_t r = 0; r < tile_height_rows; ++r)
      for (size_t c = 0; c < work_llr.cols(); ++c){
        work_llr[tile_top_row + r][c] = tile_result.tile_out[r][c];
      }
  }
}




// ===================== 顶层解码 =====================
template <typename LLR>
Matrix<LLR> ofec_decode_llr(const Matrix<LLR>& llr_mat, const Params& p,
                            std::vector<TileEarlyStopCounter>* tile_stats)
{
  const size_t N = Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;

  const size_t RROWS = llr_mat.rows();
    const size_t CCOLS = llr_mat.cols();
    if (CCOLS != N)
        throw std::invalid_argument("ofec_decode_llr: llr_mat cols != N.");

    assert(p.valid());

    const size_t TILE_HEIGHT_ROWS = p.tile_height_rows();
    const size_t TILE_STRIDE_ROWS = p.tile_stride_rows();
    const size_t WIN_HEIGHT_ROWS  = p.win_height_rows();
    const size_t POP_PUSH_ROWS    = p.pop_push_rows();
    const size_t TILES_PER_WIN    = p.TILES_PER_WIN;

    // 通道 LLR：固定（Lch）
    Matrix<LLR> channel_llr = llr_mat;

    // 工作矩阵：仅存外信息（La/Le）。首轮先验为 0。
    Matrix<LLR> work_llr(RROWS, N);
    for (size_t r = 0; r < RROWS; ++r)
        for (size_t c = 0; c < N; ++c)
            work_llr[r][c] = llr_from_float<LLR>(0.0f);

  if (RROWS < WIN_HEIGHT_ROWS) {
    if (tile_stats) {
      tile_stats->assign(p.TILES_PER_WIN, TileEarlyStopCounter{});
    }
    return channel_llr; // 不能开窗：直接回传 Lch
  }

  size_t win_start     = p.initial_win_start_rows();
  const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

  std::vector<TileEarlyStopCounter> local_tile_stats;
  std::vector<TileEarlyStopCounter>* stats_ptr = nullptr;
  if (tile_stats) {
    local_tile_stats.assign(p.TILES_PER_WIN, TileEarlyStopCounter{});
    stats_ptr = &local_tile_stats;
  }

  while (win_start <= last_ws) {
    const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;

    process_window<LLR>(work_llr, channel_llr,
                        win_start, win_end, p,
                        TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN,
                        stats_ptr);

    win_start += POP_PUSH_ROWS;
  }

  if (tile_stats) {
    *tile_stats = std::move(local_tile_stats);
  }

  // 输出：L = Lch + Le
  Matrix<LLR> out(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r)
    for (size_t c = 0; c < N; ++c) {
            const float sum = llr_to_float(channel_llr[r][c]) + llr_to_float(work_llr[r][c]);
            out[r][c] = llr_from_float<LLR>(sum);
        }
    return out;
}

// ===== 显式实例化 =====
template TileProcessResult<float>  process_tile<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool);
template TileProcessResult<int8_t> process_tile<int8_t>(const Matrix<int8_t>&,const Matrix<int8_t>&,const Params&,size_t,bool);

template void process_window<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                     size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*);
template void process_window<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&, size_t, size_t, const Params&,
                                     size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*);

template Matrix<float>  ofec_decode_llr<float >(const Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*);
template Matrix<int8_t> ofec_decode_llr<int8_t>(const Matrix<int8_t>&, const Params&, std::vector<TileEarlyStopCounter>*);
// qfloat 量化类型
template TileProcessResult<newcode::qfloat<4>> process_tile<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&,
                                                                               const Matrix<newcode::qfloat<4>>&,
                                                                               const Params&, size_t, bool);
template TileProcessResult<newcode::qfloat<5>> process_tile<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&,
                                                                               const Matrix<newcode::qfloat<5>>&,
                                                                               const Params&, size_t, bool);

template void process_window<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&, const Matrix<newcode::qfloat<4>>&, size_t, size_t, const Params&,
                                                size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*);
template void process_window<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&, const Matrix<newcode::qfloat<5>>&, size_t, size_t, const Params&,
                                                size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*);

template Matrix<newcode::qfloat<4>> ofec_decode_llr<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&, std::vector<TileEarlyStopCounter>*);
template Matrix<newcode::qfloat<5>> ofec_decode_llr<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&, std::vector<TileEarlyStopCounter>*);

} // namespace newcode
