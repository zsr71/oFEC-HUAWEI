#pragma once
#include <vector>
#include <cstdint>
#include "newcode/params.hpp"
#include "newcode/matrix.hpp"

namespace newcode {

// oFEC 编码函数：输入比特序列和系统参数，输出一个二维矩阵
Matrix<uint8_t> ofec_encode(const std::vector<uint8_t>& bits, const Params& p);

} // namespace newcode
