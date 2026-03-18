#include "new_float_only/tx/bitgen/bitgen.hpp"

#include <random>

namespace bitgen {

/**
 * 根据给定长度和随机种子生成发送端原始比特。
 * 这里不依赖解码配置，避免把单次运行参数耦合进解码核心结构。
 */
std::vector<uint8_t> generate_bits(std::size_t num_bits, int seed, bool random_bits) {
    std::vector<uint8_t> bits;
    bits.reserve(num_bits);

    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> dist(0, 1);

    for (std::size_t i = 0; i < num_bits; ++i) {
        const uint8_t bit = random_bits
                              ? static_cast<uint8_t>(dist(rng))
                              : static_cast<uint8_t>(0);
        bits.push_back(bit);
    }

    return bits;
}

} // namespace bitgen
