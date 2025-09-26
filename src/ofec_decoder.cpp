#include "newcode/ofec_decoder.hpp"
#include <stdexcept>
#include <cassert>
#include <algorithm>
#include <vector>
#include <cstddef>
#include <cstdint>

namespace newcode {

// 处理单个 Tile（占位：直通）
template <typename LLR>
Matrix<LLR> process_tile(const Matrix<LLR>& tile_in, const Params& p)
{
    (void)p;
    return tile_in;
}

// 处理一个 Window 内的所有 Tile（自下而上），内部执行 p.WINDOW_ITERS 轮
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
            for (size_t r = 0; r < tile_height_rows; ++r) {
                for (size_t c = 0; c < work_llr.cols(); ++c) {
                    tile_in[r][c] = work_llr[tile_top_row + r][c];
                }
            }

            // 2) Tile 解码（占位）
            Matrix<LLR> tile_out = process_tile<LLR>(tile_in, p);

            // 3) 回写
            for (size_t r = 0; r < tile_height_rows; ++r) {
                for (size_t c = 0; c < work_llr.cols(); ++c) {
                    work_llr[tile_top_row + r][c] = tile_out[r][c];
                }
            }
        }
    }
}

// 顶层解码：输入 LLR 矩阵，输出“更新后的 LLR 矩阵”（同类型）
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

    if (RROWS < WIN_HEIGHT_ROWS) {
        return work_llr;
    }

    size_t win_start     = p.initial_win_start_rows();
    const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

    while (win_start <= last_ws) {
        const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;

        process_window<LLR>(work_llr, win_start, win_end, p,
                            TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN);

        win_start += POP_PUSH_ROWS;
    }

    return work_llr;
}

// ===== 显式实例化（与 .hpp 中 extern template 对应）=====
// 基础类型
template Matrix<float>  process_tile<float >(const Matrix<float>&,  const Params&);
template Matrix<int8_t> process_tile<int8_t>(const Matrix<int8_t>&, const Params&);

template void process_window<float >(Matrix<float>&,  size_t, size_t, const Params&,
                                     size_t, size_t, size_t);
template void process_window<int8_t>(Matrix<int8_t>&, size_t, size_t, const Params&,
                                     size_t, size_t, size_t);

template Matrix<float>  ofec_decode_llr<float >(const Matrix<float>&,  const Params&);
template Matrix<int8_t> ofec_decode_llr<int8_t>(const Matrix<int8_t>&, const Params&);

// qfloat 量化类型（按需开启）
template Matrix<newcode::qfloat<4>> process_tile<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&);
template Matrix<newcode::qfloat<5>> process_tile<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&);

template void process_window<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&, size_t, size_t, const Params&,
                                                 size_t, size_t, size_t);
template void process_window<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&, size_t, size_t, const Params&,
                                                 size_t, size_t, size_t);

template Matrix<newcode::qfloat<4>> ofec_decode_llr<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&);
template Matrix<newcode::qfloat<5>> ofec_decode_llr<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&);

} // namespace newcode
