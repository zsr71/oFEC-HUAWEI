#include "newcode/ofec_decoder.hpp"
#include "newcode/chase256.hpp" // 保留 Chase 头；本文档内有三参前向声明
#include "newcode/bch_255_239.hpp"
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

template<typename LLR>
inline float llr_to_float(LLR x) { return static_cast<float>(x); }

template<typename LLR>
inline LLR llr_from_float(float x) {
  if constexpr (std::is_floating_point<LLR>::value) {
    return static_cast<LLR>(x);
  } else if constexpr (std::is_integral<LLR>::value && std::is_signed<LLR>::value) {
    const float lo = static_cast<float>(std::numeric_limits<LLR>::min());
    const float hi = static_cast<float>(std::numeric_limits<LLR>::max());
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return static_cast<LLR>(std::lrintf(x));
  } else {
    return LLR::from_float(x);
  }
}

template <typename LLR>
Matrix<LLR> process_tile(const Matrix<LLR>& tile_in,
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

    // 遍历底部 SBR 个 sub-block rows：从最底开始
    for (int s = 0; s < SBR; ++s)
    {
        const size_t sbr_row0_local = H - static_cast<size_t>((s + 1) * B);
        for (int r_off = 0; r_off < B; ++r_off)
        {
            const size_t row_local  = sbr_row0_local + static_cast<size_t>(r_off);
            const size_t row_global = tile_top_row_global + row_local;
            const long   R  = static_cast<long>(row_global / B);
            const int    r  = static_cast<int>(row_global % B);

            // ===== 构造两个 256 维 LLR：Lin/Lch =====
            std::array<LLR, 256> Lin256;
            std::array<LLR, 256> Lch256;

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
                    Lin256[static_cast<size_t>(k)] = llr_from_float<LLR>(Lch + La);
                    Lch256[static_cast<size_t>(k)] = llr_from_float<LLR>(Lch);
                }
                else {
                    throw std::out_of_range("process_tile: old info position out of tile range.");
                }
            }

            // (B) N..K-1，111 个"新信息"（当前行）
            for (int i = 0; i < TAKE_BITS; ++i) {
                const float Lch = llr_to_float(ch_tile[row_local][static_cast<size_t>(i)]);
                const float La  = llr_to_float(tile_in [row_local][static_cast<size_t>(i)]);
                Lin256[static_cast<size_t>(N + i)] = llr_from_float<LLR>(Lch + La);
                Lch256[static_cast<size_t>(N + i)] = llr_from_float<LLR>(Lch);
            }

            // (C) 239..254，16 个 BCH 校验位（当前行）
            for (int j = 0; j < BCH_PAR; ++j) {
                const size_t col = static_cast<size_t>(TAKE_BITS + j);
                const float Lch = llr_to_float(ch_tile[row_local][col]);
                const float La  = llr_to_float(tile_in [row_local][col]);
                Lin256[static_cast<size_t>(K + j)] = llr_from_float<LLR>(Lch + La);
                Lch256[static_cast<size_t>(K + j)] = llr_from_float<LLR>(Lch);
            }

            // (D) 255，整体偶校验位（当前行）
            {
                const size_t col = static_cast<size_t>(TAKE_BITS + BCH_PAR);
                const float Lch = llr_to_float(ch_tile[row_local][col]);
                const float La  = llr_to_float(tile_in [row_local][col]);
                Lin256[static_cast<size_t>(OVR_IDX)] = llr_from_float<LLR>(Lch + La);
                Lch256[static_cast<size_t>(OVR_IDX)] = llr_from_float<LLR>(Lch);
            }

            // ====== Chase + BCH 解码 ======
            std::array<LLR, 256> Y2{}; // 输出：Lpost - Lch
            bool produced_llr = false;

            // 可选硬判
            if (use_hard_decode)
            {
                std::array<uint8_t, 256> hard_in{};
                for (int i = 0; i < 256; ++i)
                    hard_in[static_cast<size_t>(i)] = (llr_to_float(Lin256[static_cast<size_t>(i)]) < 0.f) ? 1u : 0u;

                std::array<uint8_t, Params::BCH_N - 1> decoded{};
                if (bch_255_239_decode_hiho_cw_255(hard_in.data(), decoded.data()))
                {
                    std::array<uint8_t, Params::BCH_N> cw{};
                    for (int i = 0; i < static_cast<int>(Params::BCH_N) - 1; ++i)
                        cw[static_cast<size_t>(i)] = decoded[static_cast<size_t>(i)];
                    uint8_t parity = 0;
                    for (int i = 0; i < static_cast<int>(Params::BCH_N) - 1; ++i)
                        parity ^= cw[static_cast<size_t>(i)];
                    cw[static_cast<size_t>(OVR_IDX)] = parity;

                    const float hard_mag = std::fabs(p.HARD_LLR_MAG);
                    for (int i = 0; i < static_cast<int>(Params::BCH_N); ++i)
                    {
                        const float sign = cw[static_cast<size_t>(i)] ? -1.f : 1.f;
                        const float Lpost = sign * hard_mag;
                        const float Lch = llr_to_float(Lch256[static_cast<size_t>(i)]);
                        Y2[static_cast<size_t>(i)] = llr_from_float<LLR>(Lpost - Lch);
                    }
                    produced_llr = true;
                }
            }

            if (!use_hard_decode)
            {
                chase_decode_256<LLR>(Lin256.data(), Lch256.data(), Y2.data(), p);
                produced_llr = true;
            }

            if (!produced_llr)
                continue; // 仅选择硬判时且失败，保留原外信息
            // ====== 回写外信息到 tile_out（作为下一轮/下一组件先验） ======
            // 当前行：111 新信息 + 16 BCH + 1 overall
            for (int i = 0; i < TAKE_BITS; ++i)
                tile_out[row_local][static_cast<size_t>(i)] = Y2[static_cast<size_t>(N + i)];
            for (int j = 0; j < BCH_PAR; ++j)
                tile_out[row_local][static_cast<size_t>(TAKE_BITS + j)] = Y2[static_cast<size_t>(K + j)];
            tile_out[row_local][static_cast<size_t>(TAKE_BITS + BCH_PAR)] = Y2[static_cast<size_t>(OVR_IDX)];

            // 旧信息 128：把 Y2[0..N-1] 按映射散布回对应 (rr,cc)
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
                        Y2[static_cast<size_t>(k)];
                
            }
        } // r_off
     } // s

    return tile_out;
  } // process_tile

  // ========== 窗口处理 ==========
  template <typename LLR>
  void process_window(Matrix<LLR>& work_llr,
                      const Matrix<LLR>& channel_llr,
                      size_t win_start, size_t win_end, const Params& p,
                      size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN)
  {
    const size_t N = Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;

    for (size_t t = 0; t < TILES_PER_WIN; ++t)
    {
        const size_t tile_top_row = win_start + t * tile_stride_rows;
        const size_t tile_bottom_row = tile_top_row + tile_height_rows - 1;

        if (tile_bottom_row > win_end) break; // 越界：退出

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

        Matrix<LLR> tile_out = process_tile<LLR>(tile_in, ch_tile, tile_params,
                                                 /*tile_top_row_global=*/tile_top_row,
                                                 /*use_hard_decode=*/use_hard);

        for (size_t r = 0; r < tile_height_rows; ++r)
            for (size_t c = 0; c < work_llr.cols(); ++c){
                work_llr[tile_top_row + r][c] = tile_out[r][c];
            }
    }
  }




// ===================== 顶层解码 =====================
template <typename LLR>
Matrix<LLR> ofec_decode_llr(const Matrix<LLR>& llr_mat, const Params& p)
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
        return channel_llr; // 不能开窗：直接回传 Lch
    }

    size_t win_start     = p.initial_win_start_rows();
    const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

    while (win_start <= last_ws) {
        const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;

        process_window<LLR>(work_llr, channel_llr,win_start, win_end, p,TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN);

        win_start += POP_PUSH_ROWS;
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
template Matrix<float>  process_tile<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool);
template Matrix<int8_t> process_tile<int8_t>(const Matrix<int8_t>&,const Matrix<int8_t>&,const Params&,size_t,bool);

template void process_window<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                     size_t, size_t, size_t);
template void process_window<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&, size_t, size_t, const Params&,
                                     size_t, size_t, size_t);

template Matrix<float>  ofec_decode_llr<float >(const Matrix<float>&,  const Params&);
template Matrix<int8_t> ofec_decode_llr<int8_t>(const Matrix<int8_t>&, const Params&);
// qfloat ��������
template Matrix<newcode::qfloat<4>> process_tile<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&,
                                                                    const Matrix<newcode::qfloat<4>>&,
                                                                    const Params&, size_t, bool);
template Matrix<newcode::qfloat<5>> process_tile<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&,
                                                                    const Matrix<newcode::qfloat<5>>&,
                                                                    const Params&, size_t, bool);

template void process_window<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&, const Matrix<newcode::qfloat<4>>&, size_t, size_t, const Params&,
                                                size_t, size_t, size_t);
template void process_window<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&, const Matrix<newcode::qfloat<5>>&, size_t, size_t, const Params&,
                                                size_t, size_t, size_t);

template Matrix<newcode::qfloat<4>> ofec_decode_llr<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&);
template Matrix<newcode::qfloat<5>> ofec_decode_llr<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&);

} // namespace newcode
