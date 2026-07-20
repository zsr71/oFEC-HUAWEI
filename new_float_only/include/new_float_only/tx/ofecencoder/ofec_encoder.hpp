#pragma once
#include <vector>
#include <cstdint>
#include "new_float_only/params.hpp"
#include "new_float_only/common/matrix/matrix.hpp"

namespace ofecencoder {

// oFEC 编码函数：输入比特序列和系统参数，输出一个二维矩阵
matrix::Matrix<uint8_t> ofec_encode(const std::vector<uint8_t>& bits, const new_float_only::Params& p);

} // namespace new_float_only
