#pragma once
#include <vector>
#include <cstddef>
#include "newcode/params.hpp"
#include "newcode/ofec_llr_matrix.hpp" // Matrix<float>, Matrix<uint8_t>

namespace newcode {

// 将 LLR 矩阵按“行解码”的方式，输出一个 0/1 矩阵：
// 每行仅回填该行对应的“新信息位” (111 比特) 至列 [0 .. 110]，其余列为 0。
Matrix<uint8_t> ofec_decode_llr(const Matrix<float>& llr_mat, const Params& p);

} // namespace newcode
