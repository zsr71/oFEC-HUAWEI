#include "newcode/bitgen.hpp"
#include <random>

namespace newcode {

std::vector<uint8_t> generate_bits(const Params& params) {
    std::vector<uint8_t> bits;
    bits.reserve(params.NUM_INFO_BITS);

    std::mt19937 rng(params.BITGEN_SEED);               // 用 seed 初始化随机数
    std::uniform_int_distribution<int> dist(0, 1);

    for (size_t i = 0; i < params.NUM_INFO_BITS; i++) {
        bits.push_back(static_cast<uint8_t>(dist(rng)));
    }

    return bits;
}

} // namespace newcode
