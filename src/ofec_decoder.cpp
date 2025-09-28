#include "newcode/ofec_decoder.hpp"
#include "newcode/chase256.hpp"          // ★ 需要调用 chase_decode_256
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <cmath>
#include <type_traits>

namespace newcode {

// --------- 小工具：LLR <-> float 适配（支持 float/int8_t/qfloat<NBITS>） ---------
namespace {
template<typename LLR>
inline float llr_to_float(LLR x) { return static_cast<float>(x); }

template<typename LLR, typename std::enable_if<std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return static_cast<LLR>(x); }

template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return LLR::from_float(x); }
} // anon

// ===================== 处理单个 Tile =====================
// 实现：从 tile 底部 1/2 个 sub-block rows 提取 16/32 个码字，Chase 解码，再回写
template <typename LLR>
Matrix<LLR> process_tile(const Matrix<LLR>& tile_in, const Params& p,
                         size_t tile_top_row_global)  // ★ 新增：tile 全局顶行号
{
    const int B = static_cast<int>(Params::BITS_PER_SUBBLOCK_DIM);                 // 16
    const int N = static_cast<int>(Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM); // 128
    const int K = static_cast<int>(Params::BCH_K);                                 // 239
    const int TAKE_BITS = K - N;                                                   // 111
    const int PAR_LEN   = static_cast<int>(Params::BCH_N - Params::BCH_K);         // 16

    const size_t H = tile_in.rows();
    const size_t W = tile_in.cols(); // 应等于 N
    assert(W == static_cast<size_t>(N));

    Matrix<LLR> tile_out = tile_in; // 先拷贝，回写时覆盖相关位置

    // 需要处理的“底部 sub-block rows”个数：1 或 2（超出范围则报错）
    const int SBR = p.CHASE_SBR;
    if (SBR != 1 && SBR != 2)
        throw std::invalid_argument("process_tile: CHASE_SBR must be 1 or 2.");
    // 遍历底部 SBR 个 sub-block rows：从最底开始
    for (int s = 0; s < SBR; ++s)
    {
        const size_t sbr_row0_local = H - static_cast<size_t>((s + 1) * B); // 该 sub-block-row 的起始本地行
        // 每个 sub-block row 有 B (=16) 条 bit-rows，每条行各自对应一个码字（行内 128 + 旧信息 128-映射 + BCH/OP）
        for (int r_off = 0; r_off < B; ++r_off)
        {
            const size_t row_local  = sbr_row0_local + static_cast<size_t>(r_off);     // 本地行号
            const size_t row_global = tile_top_row_global + row_local;                 // 全局行号
            const long   R  = static_cast<long>(row_global / B);                       // 全局 block-row
            const int    r  = static_cast<int>(row_global % B);                        // 子块内行 (0..B-1)

            // ===== 构造一个 256 维码字的 LLR =====
            std::array<LLR, 256> Y256;

            // (A) 0..N-1：“旧信息”128 位（按 encoder 公式，从 tile 内其它位置取 LLR）
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
                    Y256[static_cast<size_t>(k)] =
                        tile_in[static_cast<size_t>(rr_local)][static_cast<size_t>(cc_local)];
                }
                else
                {
                    // 超出 tile 范围：报错
                    throw std::out_of_range("process_tile: old info position out of tile range.");
                }
            }

            // (B) N..K-1 = 111 个“新信息”：来自当前行的列 [0..110]
            for (int i = 0; i < TAKE_BITS; ++i)
                Y256[static_cast<size_t>(N + i)] = tile_in[row_local][static_cast<size_t>(i)];

            // (C) K..K+15：16 个 BCH parity，来自列 [111..126]
            for (int j = 0; j < PAR_LEN; ++j)
                Y256[static_cast<size_t>(K + j)] = tile_in[row_local][static_cast<size_t>(TAKE_BITS + j)];

            // (D) K+16：整体偶校验位，来自列 [127]
            Y256[static_cast<size_t>(K + PAR_LEN)] = tile_in[row_local][static_cast<size_t>(TAKE_BITS + PAR_LEN)];

            // ====== 调用 Chase–Pyndiah 取得外信息 ======
            std::array<LLR, 256> Y2;
            chase_decode_256<LLR>(Y256.data(), Y2.data(), p); // 只做流程 2~6

            // ====== 流程 7：把外信息回写到 tile ======

            // (1) 当前行的 111 个“新信息” + 16 BCH + 1 overall：直接按（B）(C)(D) 的逆映射写回
            for (int i = 0; i < TAKE_BITS; ++i)
                tile_out[row_local][static_cast<size_t>(i)] =
                    Y2[static_cast<size_t>(N + i)];

            for (int j = 0; j < PAR_LEN; ++j)
                tile_out[row_local][static_cast<size_t>(TAKE_BITS + j)] =
                    Y2[static_cast<size_t>(K + j)];

            tile_out[row_local][static_cast<size_t>(TAKE_BITS + PAR_LEN)] =
                Y2[static_cast<size_t>(K + PAR_LEN)];

            // (2) 旧信息 128：把 Y2[0..N-1] 散射回对应的 (rr,cc)
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
                    // 直接覆盖外信息（也可考虑与原 a priori 融合：如 0.5*(old + new)）
                    tile_out[static_cast<size_t>(rr_local)][static_cast<size_t>(cc_local)] =
                        Y2[static_cast<size_t>(k)];
                }
                // else: 映射到 tile 之外则跳过（由其它 tile 负责）
            }
        } // r_off
    } // s

    return tile_out;
}

// ===================== Window 内的所有 Tiles =====================
template <typename LLR>
void process_window(Matrix<LLR>& work_llr,
                    size_t win_start, size_t win_end, const Params& p,
                    size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN)
{
    assert(p.valid());

    for (size_t it = 0; it < p.WINDOW_ITERS; ++it) {
        for (size_t t = 0; t < TILES_PER_WIN; ++t) {
            const size_t tile_bottom_row = win_end - t * tile_stride_rows;
            const size_t tile_top_row    = tile_bottom_row - tile_height_rows + 1;

            // 边界防御
            assert(tile_top_row >= win_start);
            assert(tile_bottom_row < work_llr.rows());
            assert(tile_top_row + tile_height_rows - 1 <= win_end);

            // 1) 提取当前 Tile
            Matrix<LLR> tile_in(tile_height_rows, work_llr.cols());
            for (size_t r = 0; r < tile_height_rows; ++r)
                for (size_t c = 0; c < work_llr.cols(); ++c)
                    tile_in[r][c] = work_llr[tile_top_row + r][c];

            // 2) Tile 解码（使用全局顶行号传入）
            Matrix<LLR> tile_out = process_tile<LLR>(tile_in, p, /*tile_top_row_global=*/tile_top_row);
            //Matrix<LLR> tile_out = tile_in;

            // 3) 回写
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

    Matrix<LLR> work_llr = llr_mat;

    if (RROWS < WIN_HEIGHT_ROWS)
        return work_llr;

    size_t win_start     = p.initial_win_start_rows();
    const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

    while (win_start <= last_ws) {
        const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;

        process_window<LLR>(work_llr, win_start, win_end, p,TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN);

        win_start += POP_PUSH_ROWS;
    }

    return work_llr;
}

// ===== 显式实例化（与 .hpp 中 extern template 对应）=====
// 基础类型
template Matrix<float>  process_tile<float >(const Matrix<float>&,  const Params&, size_t);
template Matrix<int8_t> process_tile<int8_t>(const Matrix<int8_t>&, const Params&, size_t);

template void process_window<float >(Matrix<float>&,  size_t, size_t, const Params&,
                                     size_t, size_t, size_t);
template void process_window<int8_t>(Matrix<int8_t>&, size_t, size_t, const Params&,
                                     size_t, size_t, size_t);

template Matrix<float>  ofec_decode_llr<float >(const Matrix<float>&,  const Params&);
template Matrix<int8_t> ofec_decode_llr<int8_t>(const Matrix<int8_t>&, const Params&);

// qfloat 量化类型（按需开启）
template Matrix<newcode::qfloat<4>> process_tile<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&, size_t);
template Matrix<newcode::qfloat<5>> process_tile<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&, size_t);

template void process_window<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&, size_t, size_t, const Params&,
                                                 size_t, size_t, size_t);
template void process_window<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&, size_t, size_t, const Params&,
                                                 size_t, size_t, size_t);

template Matrix<newcode::qfloat<4>> ofec_decode_llr<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&);
template Matrix<newcode::qfloat<5>> ofec_decode_llr<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&);

} // namespace newcode
