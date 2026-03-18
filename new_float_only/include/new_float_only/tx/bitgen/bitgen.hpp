#pragma once

#include <cstddef>
#include <vector>
#include <cstdint>

namespace bitgen {

/**
 * 生成 0/1 比特流。
 * num_bits 为输出长度，seed 控制随机序列，random_bits=false 时返回全 0。
 */
std::vector<uint8_t> generate_bits(std::size_t num_bits,
                                   int seed,
                                   bool random_bits);

} // namespace bitgen
