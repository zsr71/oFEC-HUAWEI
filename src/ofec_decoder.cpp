#include "newcode/ofec_decoder.hpp"
#include "newcode/chase256.hpp" // 保留 Chase 头；本文档内有三参前向声明
#include "newcode/bch_255_239.hpp"
#include "newcode/ofec_decoder_hard.hpp"
#include "newcode/decoder_core.hpp"
#include "newcode/llr_utils.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace newcode {

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
    for (int j = 0; j < 255; ++j) {
      const float v = llr_to_float(lin_matrix[r][static_cast<size_t>(j)]);
      hard255[static_cast<size_t>(j)] = (v < 0.0f) ? 1u : 0u;
    }

    if (!bch_255_239_decode_hiho_cw_255(hard255.data(), decoded255.data()))
      return false;

    uint8_t parity255 = 0u;
    for (int j = 0; j < 255; ++j) parity255 ^= (hard255[static_cast<size_t>(j)] & 1u);
    const uint8_t overall = (llr_to_float(lin_matrix[r][255]) < 0.0f) ? 1u : 0u;

    if ((parity255 ^ overall) != 0u)
      return false;
  }

  return true;
}

namespace {

template <typename LLR>
using CoreFn = DecoderCoreResult<LLR> (*)(const Matrix<LLR>&,
                                          const Matrix<LLR>&,
                                          bool,
                                          const Params&);

template <typename LLR>
TileProcessResult<LLR> process_tile_impl(const Matrix<LLR>& tile_in,
                                         const Matrix<LLR>& ch_tile,
                                         const Params& p,
                                         size_t tile_top_row_global,
                                         bool use_hard_decode,
                                         bool normalize_extrinsic,
                                         CoreFn<LLR> core_fn)
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

  Matrix<LLR> tile_out = tile_in;

  const auto& trace_cfg = p.debug_trace;
  const bool trace_has_coords = (trace_cfg.row >= 0 && trace_cfg.col >= 0);
  const bool trace_enabled = trace_cfg.enable && trace_has_coords;
  const long trace_row = trace_has_coords ? trace_cfg.row : -1;
  const long trace_col = trace_has_coords ? trace_cfg.col : -1;
  auto should_trace = [&](bool flag, long rr, long cc) -> bool {
    return trace_enabled && flag && rr == trace_row && cc == trace_col;
  };

  const int SBR = p.CHASE_SBR;
  if (SBR != 1 && SBR != 2)
      throw std::invalid_argument("process_tile: CHASE_SBR must be 1 or 2.");

  const size_t rows_to_decode = static_cast<size_t>(SBR) * static_cast<size_t>(B);
  if (rows_to_decode == 0) {
      return TileProcessResult<LLR>{tile_out, false};
  }

  const size_t decoder_cols = static_cast<size_t>(2 * N); // 256
  Matrix<LLR> lin_matrix(rows_to_decode, decoder_cols);
  Matrix<LLR> lch_matrix(rows_to_decode, decoder_cols);
  std::vector<size_t> row_local_lookup(rows_to_decode, 0);
  std::vector<size_t> row_global_lookup(rows_to_decode, 0);

  for (int s = 0; s < SBR; ++s){
      const size_t sbr_row0_local = H - static_cast<size_t>((SBR-s) * B);
      for (int r_off = 0; r_off < B; ++r_off)
      {
          const size_t row_idx   = static_cast<size_t>(s * B + r_off);
          const size_t row_local = sbr_row0_local + static_cast<size_t>(r_off);
          const size_t row_global = tile_top_row_global + row_local;

          row_local_lookup[row_idx]  = row_local;
          row_global_lookup[row_idx] = row_global;

          const long R = static_cast<long>(row_global / static_cast<size_t>(B));
          const int  r = static_cast<int>(row_global % static_cast<size_t>(B));

          for (int k = 0; k < N; ++k)
          {
              const long br = (R ^ 1L)
                            - static_cast<long>(2 * p.NUM_GUARD_SUBROWS)
                            - static_cast<long>(2 * (N / B))
                            + static_cast<long>(2 * (k / B));
              const long bc = static_cast<long>(k / B);
              const long bit_row_in_block = static_cast<long>((k % B) ^ r);
              const long bit_col_in_block = static_cast<long>(r);

              const long rr_global = br * B + bit_row_in_block;
              const long cc_global = bc * B + bit_col_in_block;

              if (should_trace(trace_cfg.log_read_mapping, rr_global, cc_global)) {
                  std::cout << " READ Mapping k=" << k
                            << " to global pos (" << rr_global << "," << cc_global << ")" << '\n';
              }



              const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
              const long cc_local2 = cc_global;

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

          for (int i = 0; i < TAKE_BITS; ++i) {
              const int k = N + i;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t src_col = Ct * static_cast<size_t>(B) + ct;

              const float Lch = llr_to_float(ch_tile[row_local][src_col]);
              const float La  = llr_to_float(tile_in [row_local][src_col]);
              lin_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch + La);
              lch_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch);
          }

          for (int j = 0; j < BCH_PAR; ++j) {
              const int k = K + j;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t src_col = Ct * static_cast<size_t>(B) + ct;

              const float Lch = llr_to_float(ch_tile[row_local][src_col]);
              const float La  = llr_to_float(tile_in [row_local][src_col]);
              lin_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch + La);
              lch_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch);
          }

          {
              const int k = OVR_IDX;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t src_col = Ct * static_cast<size_t>(B) + ct;

              const float Lch = llr_to_float(ch_tile[row_local][src_col]);
              const float La  = llr_to_float(tile_in [row_local][src_col]);
              lin_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch + La);
              lch_matrix[row_idx][static_cast<size_t>(k)] = llr_from_float<LLR>(Lch);
          }
      }
    }

  bool early_stop_triggered = tile_should_early_stop(lin_matrix);

  auto decoder_res = core_fn(lin_matrix, lch_matrix, use_hard_decode, p);

  // Normalization only applies to soft extrinsics; skip for hard-decoded tiles.
  if (normalize_extrinsic && !use_hard_decode)
  {
      auto is_fallback = [&](float w) -> bool {
          const float target = p.beta;
          const float diff   = std::fabs(std::fabs(w) - target);
          const float tol    = 1e-4f * std::max(1.0f, target);
          return diff <= tol;
      };

      double acc = 0.0;
      std::size_t cnt = 0;

      const std::size_t Rcnt = decoder_res.lout.rows();
      const std::size_t Ccnt = decoder_res.lout.cols();

      for (std::size_t r = 0; r < Rcnt; ++r)
      {
          if (!decoder_res.produced_rows[r]) continue;
          for (std::size_t j = 0; j < Ccnt; ++j)
          {
              const float w = llr_to_float(decoder_res.lout[r][j]);
              if (is_fallback(w)) {acc += std::fabs(w/p.beta); ++cnt; continue;};
              acc += std::fabs(w);
              ++cnt;
          }
      }

      if (cnt > 0)
      {
          const float g_alpha = static_cast<float>(acc / static_cast<double>(cnt));
          if (g_alpha > 0.f)
          {
              const float scale = 1.0f / g_alpha;
              for (std::size_t r = 0; r < Rcnt; ++r)
              {
                  if (!decoder_res.produced_rows[r]) continue;
                  for (std::size_t j = 0; j < Ccnt; ++j)
                  {
                      const float w = llr_to_float(decoder_res.lout[r][j]);
                      //if (is_fallback(w)) continue;
                      decoder_res.lout[r][j] = llr_from_float<LLR>(w * scale);
                  }
              }
          }
      }
  }
  if(!use_hard_decode)
  {
      // Apply global α scaling on extrinsic outputs (moved from Chase decoder).
      const std::size_t Rcnt = decoder_res.lout.rows();
      const std::size_t Ccnt = decoder_res.lout.cols();
      for (std::size_t r = 0; r < Rcnt; ++r)
      {
          if (!decoder_res.produced_rows[r]) continue;
          for (std::size_t j = 0; j < Ccnt; ++j)
          {
              const float w = llr_to_float(decoder_res.lout[r][j]);
              decoder_res.lout[r][j] = llr_from_float<LLR>(w * p.ALPHA);
          }
      }
  }

  for (int s = 0; s < SBR; ++s)
  {
      for (int r_off = 0; r_off < B; ++r_off)
      {
          const size_t row_idx = static_cast<size_t>(s * B + r_off);
          if (row_idx >= decoder_res.produced_rows.size()) continue;
          if (!decoder_res.produced_rows[row_idx]) continue;

          const size_t row_local  = row_local_lookup[row_idx];
          const size_t row_global = row_global_lookup[row_idx];
          const auto&  lout_row   = decoder_res.lout[row_idx];

          const int r = static_cast<int>(row_global % static_cast<size_t>(B));

          for (int i = 0; i < TAKE_BITS; ++i) {
              const int k = N + i;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t col = Ct * static_cast<size_t>(B) + ct;
              tile_out[row_local][col] = lout_row[static_cast<size_t>(k)];
          }
          for (int j = 0; j < BCH_PAR; ++j) {
              const int k = K + j;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t col = Ct * static_cast<size_t>(B) + ct;
              tile_out[row_local][col] = lout_row[static_cast<size_t>(k)];
          }
          {
              const int k = OVR_IDX;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t col = Ct * static_cast<size_t>(B) + ct;
              tile_out[row_local][col] = lout_row[static_cast<size_t>(k)];
          }

          const long R = static_cast<long>(row_global / static_cast<size_t>(B));
          for (int k = 0; k < N; ++k)
          {
              const long br = (R ^ 1L)
                            - static_cast<long>(2 * p.NUM_GUARD_SUBROWS)
                            - static_cast<long>(2 * (N / B))
                            + static_cast<long>(2 * (k / B));
              const long bc = static_cast<long>(k / B);
              const long bit_row_in_block = static_cast<long>((k % B) ^ r);
              const long bit_col_in_block = static_cast<long>(r);

              const long rr_global = br * B + bit_row_in_block;
              const long cc_global = bc * B + bit_col_in_block;

              if (should_trace(trace_cfg.log_write_mapping, rr_global, cc_global)) {
                  std::cout << " WRITE Mapping k=" << k
                            << " to global pos (" << rr_global << "," << cc_global << ")" << '\n';
                  std::cout << "  Value=" << llr_to_float(lout_row[static_cast<size_t>(k)]) << '\n';
              }


              const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
              const long cc_local2 = cc_global;

              const bool in_range =
                  (rr_local2 >= 0 && rr_local2 < static_cast<long>(H) &&
                   cc_local2 >= 0 && cc_local2 < static_cast<long>(W));

              assert(in_range && "process_tile: write-back out of tile range");

              tile_out[static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)] =
                  lout_row[static_cast<size_t>(k)];
          }
      }
  }

  return TileProcessResult<LLR>{std::move(tile_out), early_stop_triggered};
}

template <typename LLR>
void process_window_impl(Matrix<LLR>& work_llr,
                         const Matrix<LLR>& channel_llr,
                         size_t win_start, size_t win_end, const Params& p,
                         size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                         std::vector<TileEarlyStopCounter>* tile_stats,
                         bool normalize_extrinsic,
                         CoreFn<LLR> core_fn)
{
  (void)win_start;
  const auto& trace_cfg = p.debug_trace;
  const bool trace_has_coords = (trace_cfg.row >= 0 && trace_cfg.col >= 0);
  const bool trace_mismatch =
      trace_cfg.enable && trace_cfg.log_mismatch && trace_has_coords;
  const long trace_row = trace_has_coords ? trace_cfg.row : -1;
  const long trace_col = trace_has_coords ? trace_cfg.col : -1;
  const size_t N = Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;

  for (size_t t = 0; t < TILES_PER_WIN; ++t)
  {
        const size_t tile_bottom_row = win_end  - t * tile_stride_rows;
        const size_t tile_top_row    = tile_bottom_row + 1 - tile_height_rows;

        const size_t tile_height_rows_actual = tile_bottom_row - tile_top_row + 1;
        Matrix<LLR> tile_in(tile_height_rows_actual, N);
        Matrix<LLR> ch_tile(tile_height_rows_actual, N);

        for (size_t r = 0; r < tile_height_rows_actual; ++r) {
            for (size_t c = 0; c < N; ++c) {
                tile_in[r][c]  = work_llr[tile_top_row + r][c];
                ch_tile[r][c]  = channel_llr[tile_top_row + r][c];
            }
        }

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

    TileProcessResult<LLR> tile_result = process_tile_impl<LLR>(tile_in, ch_tile, tile_params,
                                                                /*tile_top_row_global=*/tile_top_row,
                                                                /*use_hard_decode=*/use_hard,
                                                                /*normalize_extrinsic=*/normalize_extrinsic,
                                                                core_fn);

    if (tile_stats && t < tile_stats->size()) {
      auto& counter = (*tile_stats)[t];
      counter.total += 1;
      if (tile_result.early_stop_triggered) {
        counter.triggered += 1;
      }
    }

    for (size_t r = 0; r < tile_height_rows_actual; ++r) {
      const size_t global_row = tile_top_row + r;
      for (size_t c = 0; c < work_llr.cols(); ++c) {
        const auto incoming = tile_result.tile_out[r][c];
        if (trace_mismatch &&
            static_cast<long>(global_row) == trace_row &&
            static_cast<long>(c) == trace_col) {
          const float existing_val = llr_to_float(work_llr[global_row][c]);
          const float incoming_val = llr_to_float(incoming);
          if (existing_val != incoming_val) {
            std::cout << "Mismatch at work_llr[" << trace_row << "][" << trace_col
                      << "]: tile index " << t
                      << " incoming=" << incoming_val << '\n';
          }
        }
        work_llr[global_row][c] = incoming;
      }
    }
  }
}

template <typename LLR>
Matrix<LLR> ofec_decode_llr_impl(const Matrix<LLR>& llr_mat, const Params& p,
                                 std::vector<TileEarlyStopCounter>* tile_stats,
                                 bool normalize_extrinsic,
                                 CoreFn<LLR> core_fn)
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

  Matrix<LLR> channel_llr = llr_mat;

  Matrix<LLR> work_llr(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
          work_llr[r][c] = llr_from_float<LLR>(0.0f);

  if (RROWS < WIN_HEIGHT_ROWS) {
    if (tile_stats) {
      tile_stats->assign(p.TILES_PER_WIN, TileEarlyStopCounter{});
    }
    return channel_llr;
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

    process_window_impl<LLR>(work_llr, channel_llr,
                             win_start, win_end, p,
                             TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN,
                             stats_ptr,
                             normalize_extrinsic,
                             core_fn);

    win_start += POP_PUSH_ROWS;
  }

  if (tile_stats) {
    *tile_stats = std::move(local_tile_stats);
  }

  Matrix<LLR> out(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r) {
    for (size_t c = 0; c < N; ++c) {
      const float sum = llr_to_float(channel_llr[r][c]) + llr_to_float(work_llr[r][c]);
      out[r][c] = llr_from_float<LLR>(sum);
    }
  }
  return out;
}

} // namespace

template <typename LLR>
TileProcessResult<LLR> process_tile_plain(const Matrix<LLR>& tile_in,
                                          const Matrix<LLR>& ch_tile,
                                          const Params& p,
                                          size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic)
{
  return process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                           use_hard_decode, normalize_extrinsic,
                           &Decoder_Core_plain<LLR>);
}

template <typename LLR>
TileProcessResult<LLR> process_tile_ebchPF(const Matrix<LLR>& tile_in,
                                           const Matrix<LLR>& ch_tile,
                                           const Params& p,
                                           size_t tile_top_row_global,
                                           bool use_hard_decode,
                                           bool normalize_extrinsic)
{
  return process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                           use_hard_decode, normalize_extrinsic,
                           &Decoder_Core_ebchPF<LLR>);
}

template <typename LLR>
void process_window_plain(Matrix<LLR>& work_llr,
                          const Matrix<LLR>& channel_llr,
                          size_t win_start, size_t win_end, const Params& p,
                          size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                          std::vector<TileEarlyStopCounter>* tile_stats,
                          bool normalize_extrinsic)
{
  process_window_impl(work_llr, channel_llr, win_start, win_end, p,
                      tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                      tile_stats, normalize_extrinsic,
                      &Decoder_Core_plain<LLR>);
}

template <typename LLR>
void process_window_ebchPF(Matrix<LLR>& work_llr,
                           const Matrix<LLR>& channel_llr,
                           size_t win_start, size_t win_end, const Params& p,
                           size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                           std::vector<TileEarlyStopCounter>* tile_stats,
                           bool normalize_extrinsic)
{
  process_window_impl(work_llr, channel_llr, win_start, win_end, p,
                      tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                      tile_stats, normalize_extrinsic,
                      &Decoder_Core_ebchPF<LLR>);
}

template <typename LLR>
Matrix<LLR> ofec_decode_llr_plain(const Matrix<LLR>& llr_mat, const Params& p,
                                  std::vector<TileEarlyStopCounter>* tile_stats,
                                  bool normalize_extrinsic)
{
  return ofec_decode_llr_impl(llr_mat, p, tile_stats, normalize_extrinsic,
                              &Decoder_Core_plain<LLR>);
}

template <typename LLR>
Matrix<LLR> ofec_decode_llr_ebchPF(const Matrix<LLR>& llr_mat, const Params& p,
                                   std::vector<TileEarlyStopCounter>* tile_stats,
                                   bool normalize_extrinsic)
{
  return ofec_decode_llr_impl(llr_mat, p, tile_stats, normalize_extrinsic,
                              &Decoder_Core_ebchPF<LLR>);
}

// ===== 显式实例化 =====
template TileProcessResult<float>  process_tile_plain<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool, bool);
template TileProcessResult<int8_t> process_tile_plain<int8_t>(const Matrix<int8_t>&, const Matrix<int8_t>&, const Params&, size_t, bool, bool);
template TileProcessResult<float>  process_tile_ebchPF<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool, bool);
template TileProcessResult<int8_t> process_tile_ebchPF<int8_t>(const Matrix<int8_t>&, const Matrix<int8_t>&, const Params&, size_t, bool, bool);

template void process_window_plain<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool);
template void process_window_plain<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&, size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool);
template void process_window_ebchPF<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool);
template void process_window_ebchPF<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&, size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool);

template Matrix<float>  ofec_decode_llr_plain<float >(const Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool);
template Matrix<int8_t> ofec_decode_llr_plain<int8_t>(const Matrix<int8_t>&, const Params&, std::vector<TileEarlyStopCounter>*, bool);
template Matrix<float>  ofec_decode_llr_ebchPF<float >(const Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool);
template Matrix<int8_t> ofec_decode_llr_ebchPF<int8_t>(const Matrix<int8_t>&, const Params&, std::vector<TileEarlyStopCounter>*, bool);

#define INSTANTIATE_QFLOAT(N) \
template TileProcessResult<newcode::qfloat<N>> process_tile_plain<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    const Params&, size_t, bool, bool); \
template TileProcessResult<newcode::qfloat<N>> process_tile_ebchPF<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    const Params&, size_t, bool, bool); \
template void process_window_plain<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    size_t, size_t, const Params&, \
    size_t, size_t, size_t, \
    std::vector<TileEarlyStopCounter>*, bool); \
template void process_window_ebchPF<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    size_t, size_t, const Params&, \
    size_t, size_t, size_t, \
    std::vector<TileEarlyStopCounter>*, bool); \
template Matrix<newcode::qfloat<N>> ofec_decode_llr_plain<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Params&, \
    std::vector<TileEarlyStopCounter>*, bool); \
template Matrix<newcode::qfloat<N>> ofec_decode_llr_ebchPF<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Params&, \
    std::vector<TileEarlyStopCounter>*, bool); \
template void process_window<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    std::size_t, std::size_t, const Params&, \
    std::size_t, std::size_t, std::size_t, \
    std::vector<TileEarlyStopCounter>*, bool);

INSTANTIATE_QFLOAT(2)
INSTANTIATE_QFLOAT(3)
INSTANTIATE_QFLOAT(4)
INSTANTIATE_QFLOAT(5)
INSTANTIATE_QFLOAT(6)
INSTANTIATE_QFLOAT(7)
INSTANTIATE_QFLOAT(8)
INSTANTIATE_QFLOAT(9)
INSTANTIATE_QFLOAT(10)
INSTANTIATE_QFLOAT(11)
INSTANTIATE_QFLOAT(12)
INSTANTIATE_QFLOAT(13)
INSTANTIATE_QFLOAT(14)
INSTANTIATE_QFLOAT(15)

#undef INSTANTIATE_QFLOAT

template void process_window<float >(Matrix<float>&,  const Matrix<float>&,
                                     std::size_t, std::size_t, const Params&,
                                     std::size_t, std::size_t, std::size_t,
                                     std::vector<TileEarlyStopCounter>*, bool);
template void process_window<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&,
                                     std::size_t, std::size_t, const Params&,
                                     std::size_t, std::size_t, std::size_t,
                                     std::vector<TileEarlyStopCounter>*, bool);

} // namespace newcode
