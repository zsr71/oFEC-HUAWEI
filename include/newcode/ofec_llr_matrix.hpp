#pragma once
#include <vector>
#include <stdexcept>

// 如果 Matrix 的定义不在这个头里，请按你的工程实际修改路径
#include "newcode/matrix.hpp"

namespace newcode {

/**
 * @brief 将一维 LLR（row-major）还原为 rows x cols 的 Matrix<float>
 * @throws std::invalid_argument 当 llr.size() != rows*cols 时
 */
Matrix<float> llr_to_matrix_row_major(const std::vector<float>& llr,
                                      size_t rows,
                                      size_t cols);

/**
 * @brief 将一维 LLR（row-major）填充到已分配的 Matrix<float>（尺寸需匹配）
 * @throws std::invalid_argument 当 llr.size() != M.rows()*M.cols() 时
 */
void fill_llr_matrix_row_major(const std::vector<float>& llr,
                               Matrix<float>& M);

} // namespace newcode
