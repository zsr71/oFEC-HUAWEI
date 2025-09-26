#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include "newcode/matrix.hpp"
#include "newcode/params.hpp"

namespace newcode {

// ====== 顶层解码：输入 LLR 矩阵，输出“更新后的 LLR 矩阵” ======
// 支持任意 LLR 标量类型（float、int8_t 等），用于软/定点流程复用。
template <typename LLR>
Matrix<LLR> ofec_decode_llr(const Matrix<LLR>& llr_mat, const Params& p);

// ====== 窗口处理：对 work_llr 的 [win_start, win_end] 做一轮（或多轮）tile 扫描 ======
// 与当前 cpp 实现保持一致：就地修改，不返回副本。
template <typename LLR>
void process_window(Matrix<LLR>& work_llr,
                    std::size_t win_start, std::size_t win_end, const Params& p,
                    std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN);

// ====== Tile 处理：对单个 tile 的 LLR 做一次“占位/解码/更新”（目前可直通） ======
template <typename LLR>
Matrix<LLR> process_tile(const Matrix<LLR>& tile_in, const Params& p);

// ——（可选）为常用类型做显式外部声明，减少编译时间 ——
// 这些会在 ofec_decoder.cpp 里用 `template ...` 显式实例化。
extern template Matrix<float>  ofec_decode_llr<float >(const Matrix<float>&,  const Params&);
extern template Matrix<int8_t> ofec_decode_llr<int8_t>(const Matrix<int8_t>&, const Params&);

extern template void process_window<float >(Matrix<float>&,  std::size_t, std::size_t, const Params&,
                                            std::size_t, std::size_t, std::size_t);
extern template void process_window<int8_t>(Matrix<int8_t>&, std::size_t, std::size_t, const Params&,
                                            std::size_t, std::size_t, std::size_t);

extern template Matrix<float>  process_tile<float >(const Matrix<float>&,  const Params&);
extern template Matrix<int8_t> process_tile<int8_t>(const Matrix<int8_t>&, const Params&);

} // namespace newcode
