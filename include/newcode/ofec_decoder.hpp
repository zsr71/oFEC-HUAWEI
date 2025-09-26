#pragma once
#include <cstddef>
#include <cstdint>
#include <vector>
#include "newcode/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"   // 为 extern template 可见类型声明

namespace newcode {

// 顶层解码：输入 LLR 矩阵，输出“更新后的 LLR 矩阵”（同类型）
template <typename LLR>
Matrix<LLR> ofec_decode_llr(const Matrix<LLR>& llr_mat, const Params& p);

// 窗口处理：对 work_llr 的 [win_start, win_end] 做一轮（或多轮）tile 扫描（就地修改）
template <typename LLR>
void process_window(Matrix<LLR>& work_llr,
                    std::size_t win_start, std::size_t win_end, const Params& p,
                    std::size_t tile_height_rows, std::size_t tile_stride_rows, std::size_t TILES_PER_WIN);

// Tile 处理：对单个 tile 的 LLR 做一次“占位/解码/更新”
template <typename LLR>
Matrix<LLR> process_tile(const Matrix<LLR>& tile_in, const Params& p);

// ===== extern template（减少重复实例化）=====
// 基础类型
extern template Matrix<float>  ofec_decode_llr<float >(const Matrix<float>&,  const Params&);
extern template Matrix<int8_t> ofec_decode_llr<int8_t>(const Matrix<int8_t>&, const Params&);
extern template void process_window<float >(Matrix<float>&,  std::size_t, std::size_t, const Params&,
                                            std::size_t, std::size_t, std::size_t);
extern template void process_window<int8_t>(Matrix<int8_t>&, std::size_t, std::size_t, const Params&,
                                            std::size_t, std::size_t, std::size_t);
extern template Matrix<float>  process_tile<float >(const Matrix<float>&,  const Params&);
extern template Matrix<int8_t> process_tile<int8_t>(const Matrix<int8_t>&, const Params&);

// qfloat 量化类型（按需开启）
extern template Matrix<newcode::qfloat<4>> ofec_decode_llr<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&);
extern template Matrix<newcode::qfloat<5>> ofec_decode_llr<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&);
extern template void process_window<newcode::qfloat<4>>(Matrix<newcode::qfloat<4>>&, std::size_t, std::size_t, const Params&,
                                                        std::size_t, std::size_t, std::size_t);
extern template void process_window<newcode::qfloat<5>>(Matrix<newcode::qfloat<5>>&, std::size_t, std::size_t, const Params&,
                                                        std::size_t, std::size_t, std::size_t);
extern template Matrix<newcode::qfloat<4>> process_tile<newcode::qfloat<4>>(const Matrix<newcode::qfloat<4>>&, const Params&);
extern template Matrix<newcode::qfloat<5>> process_tile<newcode::qfloat<5>>(const Matrix<newcode::qfloat<5>>&, const Params&);

} // namespace newcode
