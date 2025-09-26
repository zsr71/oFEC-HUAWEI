#pragma once
#include <cstddef>
#include <cstdint>

namespace newcode {

template<typename T> class Matrix;
struct Params;

// 将浮点 LLR 矩阵量化为有符号整型（存储位宽由 p.LLR_BITS 决定，常用 4/5）。
// 结果范围为 [-Q, +Q]，Q = 2^(LLR_BITS-1)-1，存放在 int8_t 中。
Matrix<int8_t> quantize_llr_to_int8(const Matrix<float>& in, const Params& p);

// （可选）反量化：把量化的整型 LLR 还原成 float（便于可视化或复用浮点接口）
Matrix<float> dequantize_llr_to_float(const Matrix<int8_t>& in, const Params& p);

// （工具）仅把 int8_t LLR 做类型转换为 float（不做缩放/反量化，符号阈值保持一致）
Matrix<float> cast_qllr_to_float(const Matrix<int8_t>& in);

} // namespace newcode
