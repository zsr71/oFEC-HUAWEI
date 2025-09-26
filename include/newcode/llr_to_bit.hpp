#pragma once
#include <cstdint>
#include <cstddef>

namespace newcode {

// 前向声明：避免头文件强依赖 Matrix 的定义文件
template<typename T> class Matrix;

// 将 LLR 矩阵转换为比特矩阵（LLR>=0 -> 0, 否则 1）
Matrix<uint8_t> llr_to_bit(const Matrix<float>& llr_mat);

} // namespace newcode
