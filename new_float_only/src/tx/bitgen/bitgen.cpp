#include "new_float_only/tx/bitgen/bitgen.hpp"

#include <random>

namespace bitgen {

/**
 * 根据给定长度和随机种子生成发送端原始比特。
 * 这里不依赖解码配置，避免把单次运行参数耦合进解码核心结构。
 */
std::vector<uint8_t> generate_bits(std::size_t num_bits, int seed, bool random_bits) {
    // 预留足够容量，避免逐个 push 时频繁扩容。
    std::vector<uint8_t> bits;
    bits.reserve(num_bits);

    // 固定种子的伪随机数引擎，保证同一 seed 下结果可复现。
    std::mt19937 rng(seed);
    // 发送端比特只有 0/1 两种取值，因此直接生成闭区间 [0, 1]。
    std::uniform_int_distribution<int> dist(0, 1);

    for (std::size_t i = 0; i < num_bits; ++i) {
        // `random_bits=false` 时生成全 0 比特流，便于做对照实验或定向排查。
        const uint8_t bit = random_bits
                              ? static_cast<uint8_t>(dist(rng))
                              : static_cast<uint8_t>(0);
        bits.push_back(bit);
    }

    return bits;
}

} // namespace bitgen
