#pragma once

#include <cstddef>
#include "newcode/matrix.hpp"
#include "newcode/params.hpp"

namespace newcode {

// 将前缀部分（已知为 0 的保护/信息子行）覆盖为强置信度的 “bit=0” LLR。
// 这在接收端可利用先验知识，避免解码器在这些位置反复迭代。
void apply_known_zero_prefix(Matrix<float>& llr_mat, const Params& p);

} // namespace newcode

