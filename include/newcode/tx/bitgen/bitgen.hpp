#pragma once
#include <vector>
#include <cstdint>
#include "newcode/params.hpp"

namespace bitgen {

// 根据 Params 生成随机比特（0 或 1）
std::vector<uint8_t> generate_bits(const newcode::Params& params);

} // namespace bitgen
