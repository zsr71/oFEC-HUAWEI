#pragma once
#include <vector>
#include <cstdint>
#include "newcode/params.hpp"
#include "newcode/common/matrix/matrix.hpp"

namespace ofecencoder {

// oFEC 编码函数：输入比特序列和系统参数，输出一个二维矩阵
matrix::Matrix<uint8_t> ofec_encode(const std::vector<uint8_t>& bits, const newcode::Params& p);

} // namespace newcode
