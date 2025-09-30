#include "newcode/ofec_decoder.hpp"
#include "newcode/chase256.hpp"          // chase_decode_256: 返回纯外信息 L_e
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <type_traits>

namespace newcode {

namespace {
  template<typename LLR>
  inline float llr_to_float(LLR x) { return static_cast<float>(x); }

  template<typename LLR, typename std::enable_if<std::is_arithmetic<LLR>::value, int>::type = 0>
  inline LLR llr_from_float(float x) { return static_cast<LLR>(x); }

  template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
  inline LLR llr_from_float(float x) { return LLR::from_float(x); }
}

// ===================== 处理单个 Tile =====================
template <typename LLR>
Matrix<LLR> process_tile(const Matrix<LLR>& tile_in,
                         const Matrix<LLR>& ch_tile,
                         const Params& p,
                         size_t tile_top_row_global)
{
    const int B = static_cast<int>(Params::BITS_PER_SUBBLOCK_DIM);                 // 16
    const int N = static_cast<int>(Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM); // 128
    const int K = static_cast<int>(Params::BCH_K);                                 // 239
    const int TAKE_BITS = K - N;                                                   // 111
    const int PAR_LEN   = static_cast<int>(Params::BCH_N - Params::BCH_K);         // 16

    const size_t H = tile_in.rows();
    const size_t W = tile_in.cols(); // 应等于 N
    assert(W == static_cast<size_t>(N));
    assert(ch_tile.rows() == H && ch_tile.cols() == W);

    Matrix<LLR> tile_out = tile_in; // 先拷贝，回写外信息时覆盖相关位置

    // 需要处理的“底部 sub-block rows”个数：1 或 2（超出范围则报错）
    const int SBR = p.CHASE_SBR;
    if (SBR != 1 && SBR != 2)
        throw std::invalid_argument("process_tile: CHASE_SBR must be 1 or 2.");

    // 遍历底部 SBR 个 sub-block rows：从最底开始
    for (int s = 0; s < SBR; ++s)
    {
        const size_t sbr_row0_local = H - static_cast<size_t>((s + 1) * B); // 该 sub-block-row 的起始本地行
        // 每个 sub-block row 有 B (=16) 条 bit-rows
        for (int r_off = 0; r_off < B; ++r_off)
        {
            const size_t row_local  = sbr_row0_local + static_cast<size_t>(r_off);     // 本地行号
            const size_t row_global = tile_top_row_global + row_local;                 // 全局行号
            const long   R  = static_cast<long>(row_global / B);                       // 全局 block-row
            const int    r  = static_cast<int>(row_global % B);                        // 子块内行 (0..B-1)

            // ===== 构造一个 256 维码字的“输入” LLR：Y256 = L_ch + L_a =====
            std::array<LLR, 256> Y256;

            // (A) 0..N-1：“旧信息”128 位
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

                // 换算到本 tile 的本地行列
                const long rr_local = rr_global - static_cast<long>(tile_top_row_global);
                const long cc_local = cc_global; // 列方向本来就是 [0..N-1]

                if (rr_local >= 0 && rr_local < static_cast<long>(H) &&
                    cc_local >= 0 && cc_local < static_cast<long>(W))
                {
                    const float Yin = llr_to_float(ch_tile[static_cast<size_t>(rr_local)][static_cast<size_t>(cc_local)])
                                     + llr_to_float(tile_in [static_cast<size_t>(rr_local)][static_cast<size_t>(cc_local)]);
                    Y256[static_cast<size_t>(k)] = llr_from_float<LLR>(Yin);
                }
                else
                {
                    throw std::out_of_range("process_tile: old info position out of tile range.");
                }
            }

            // (B) N..K-1 = 111 个“新信息”
            for (int i = 0; i < TAKE_BITS; ++i) {
                const float Yin = llr_to_float(ch_tile[row_local][static_cast<size_t>(i)])
                                 + llr_to_float(tile_in [row_local][static_cast<size_t>(i)]);
                Y256[static_cast<size_t>(N + i)] = llr_from_float<LLR>(Yin);
            }

            // (C) K..K+15：16 个 BCH parity
            for (int j = 0; j < PAR_LEN; ++j) {
                const size_t cc = static_cast<size_t>(TAKE_BITS + j);
                const float Yin = llr_to_float(ch_tile[row_local][cc])
                                 + llr_to_float(tile_in [row_local][cc]);
                Y256[static_cast<size_t>(K + j)] = llr_from_float<LLR>(Yin);
            }

            // (D) K+16：整体偶校验位
            {
                const size_t cc = static_cast<size_t>(TAKE_BITS + PAR_LEN);
                const float Yin = llr_to_float(ch_tile[row_local][cc])
                                 + llr_to_float(tile_in [row_local][cc]);
                Y256[static_cast<size_t>(K + PAR_LEN)] = llr_from_float<LLR>(Yin);
            }

            // ====== 调用 Chase–Pyndiah：返回“纯外信息” Y2 = L_e ======
            std::array<LLR, 256> Y2{}; // ★ 改为零初始化
            chase_decode_256<LLR>(Y256.data(), Y2.data(), p); // 输出 L_e

            // ====== 写回外信息到 tile_out（作为下一轮/下一组件的先验 L_a） ======
            // (1) 当前行的 111 个“新信息” + 16 BCH + 1 overall
            for (int i = 0; i < TAKE_BITS; ++i)
                tile_out[row_local][static_cast<size_t>(i)] = Y2[static_cast<size_t>(N + i)];
            for (int j = 0; j < PAR_LEN; ++j)
                tile_out[row_local][static_cast<size_t>(TAKE_BITS + j)] = Y2[static_cast<size_t>(K + j)];
            tile_out[row_local][static_cast<size_t>(TAKE_BITS + PAR_LEN)] = Y2[static_cast<size_t>(K + PAR_LEN)];

            // (2) 旧信息 128：把 Y2[0..N-1] 按映射散射回对应 (rr,cc)
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

                const long rr_local = rr_global - static_cast<long>(tile_top_row_global);
                const long cc_local = cc_global;

                if (rr_local >= 0 && rr_local < static_cast<long>(H) &&
                    cc_local >= 0 && cc_local < static_cast<long>(W))
                {
                    tile_out[static_cast<size_t>(rr_local)][static_cast<size_t>(cc_local)] =
                        Y2[static_cast<size_t>(k)];
                }
            }
        } // r_off
    } // s

    return tile_out; // 外信息矩阵切片（L_e）
}

// ===================== Window 内的所有 Tiles =====================
template <typename LLR>
void process_window(Matrix<LLR>& work_llr,
                    const Matrix<LLR>& channel_llr,
                    size_t win_start, size_t win_end, const Params& p,
                    size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN)
{
    assert(p.valid());

    for (size_t it = 0; it < p.WINDOW_ITERS; ++it) {
        Params p_it = p;
        p_it.CP_ITER = static_cast<int>(it); // ★ 把迭代号写入，供 chase256 选系数

        for (size_t t = 0; t < TILES_PER_WIN; ++t) {
            const size_t tile_bottom_row = win_end - t * tile_stride_rows;
            const size_t tile_top_row    = tile_bottom_row - tile_height_rows + 1;

            // 边界防御
            assert(tile_top_row >= win_start);
            assert(tile_bottom_row < work_llr.rows());
            assert(tile_top_row + tile_height_rows - 1 <= win_end);

            // 1) 提取当前 Tile：先验/外信息 与 通道 LLR 两个切片
            Matrix<LLR> tile_in(tile_height_rows, work_llr.cols());
            Matrix<LLR> ch_tile (tile_height_rows, work_llr.cols());
            for (size_t r = 0; r < tile_height_rows; ++r) {
                for (size_t c = 0; c < work_llr.cols(); ++c) {
                    tile_in[r][c] = work_llr    [tile_top_row + r][c];
                    ch_tile [r][c] = channel_llr[tile_top_row + r][c];
                }
            }

            // 2) Tile 解码（传入通道切片 + 先验切片）
            Matrix<LLR> tile_out = process_tile<LLR>(tile_in, ch_tile, p_it, /*tile_top_row_global=*/tile_top_row);

            // 3) 回写外信息
            for (size_t r = 0; r < tile_height_rows; ++r)
                for (size_t c = 0; c < work_llr.cols(); ++c)
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

    // 通道 LLR：固定不变（L_ch）
    Matrix<LLR> channel_llr = llr_mat;

    // 工作矩阵：仅存外信息（L_a/L_e）。首轮先验为 0。
    Matrix<LLR> work_llr(RROWS, N);
    for (size_t r = 0; r < RROWS; ++r)
        for (size_t c = 0; c < N; ++c)
            work_llr[r][c] = llr_from_float<LLR>(0.0f);

    if (RROWS < WIN_HEIGHT_ROWS) {
        // 没法开窗迭代：直接返回 L_ch（等价 L_ch + 0）
        return channel_llr;
    }

    size_t win_start     = p.initial_win_start_rows();
    const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

    while (win_start <= last_ws) {
        const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;

        process_window<LLR>(work_llr, channel_llr,
                            win_start, win_end, p,
                            TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN);

        win_start += POP_PUSH_ROWS;
    }

    // === 输出：L_ch + L_e ===
    Matrix<LLR> out(RROWS, N);
    for (size_t r = 0; r < RROWS; ++r)
        for (size_t c = 0; c < N; ++c) {
            const float sum = llr_to_float(channel_llr[r][c]) + llr_to_float(work_llr[r][c]);
            out[r][c] = llr_from_float<LLR>(sum);
        }
    return out;
}

// ===== 显式实例化 =====
template Matrix<float>  process_tile<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t);
template Matrix<int8_t> process_tile<int8_t>(const Matrix<int8_t>&, const Matrix<int8_t>&, const Params&, size_t);

template void process_window<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                     size_t, size_t, size_t);
template void process_window<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&, size_t, size_t, const Params&,
                                     size_t, size_t, size_t);

template Matrix<float>  ofec_decode_llr<float >(const Matrix<float>&,  const Params&);
template Matrix<int8_t> ofec_decode_llr<int8_t>(const Matrix<int8_t>&, const Params&);

// qfloat 量化类型
template Matrix<newcode::qfloat<4>> process_tile<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&,
                                                                     const Matrix<newcode::qfloat<4>>&,
                                                                     const Params&, size_t);
template Matrix<newcode::qfloat<5>> process_tile<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&,
                                                                     const Matrix<newcode::qfloat<5>>&,
                                                                     const Params&, size_t);

template void process_window<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&, const Matrix<newcode::qfloat<4>>&, size_t, size_t, const Params&,
                                                 size_t, size_t, size_t);
template void process_window<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&, const Matrix<newcode::qfloat<5>>&, size_t, size_t, const Params&,
                                                 size_t, size_t, size_t);

template Matrix<newcode::qfloat<4>> ofec_decode_llr<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&);
template Matrix<newcode::qfloat<5>> ofec_decode_llr<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&);

} // namespace newcode
