#include "new_float_only/common/bch/bch_255_239.hpp"
#include <algorithm>

namespace bch
{
// ---------------- LFSR 奇偶（保留你的实现） ----------------
static inline std::array<uint8_t,16> parity_core_239(const uint8_t* info239)
{
    std::array<uint8_t,16> reg{}; // 全 0
    for (int i = 239 - 1; i >= 0; --i)
    {
        const uint8_t feedback = (info239[i] & 1u) ^ reg[15];
        for (int j = 15; j > 0; --j)
            reg[j] = static_cast<uint8_t>( reg[j - 1] ^ (bch::G_COEFFS[j] & feedback) );
        reg[0] = static_cast<uint8_t>( bch::G_COEFFS[0] & feedback ); // g0=1 => reg[0]=feedback
    }
    return reg; // 16 位校验
}

std::array<uint8_t,16> bch_255_239_parity(const std::vector<uint8_t>& info239)
{
    std::array<uint8_t,239> buf{};
    const int upto = static_cast<int>(std::min<size_t>(239, info239.size()));
    for (int i = 0; i < upto; ++i) buf[i] = (info239[i] & 1u);
    return parity_core_239(buf.data());
}

} // namespace bch
