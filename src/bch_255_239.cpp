#include "newcode/bch_255_239.hpp"
#include <array>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cstddef>

namespace newcode
{
// 生成多项式： g(y) = y^16 + y^14 + y^13 + y^11 + y^10 + y^9 + y^8 + y^6 + y^5 + y + 1
// LFSR 长度 m=16，这里只存 0..15 次项（最高次 y^16 隐含为 1）
static constexpr uint8_t G_COEFFS[16] = {
    /*y^0..y^15*/
    1, 1, 0, 0, 0, 1, 1, 0,
    1, 1, 1, 1, 0, 1, 1, 0
};

// 内部：基于 LFSR 的多项式除法，输入 239 位（MSB 优先：i=238..0）
static inline std::array<uint8_t,16> parity_core_239(const uint8_t* info239)
{
    std::array<uint8_t,16> reg{}; // 全 0
    for (int i = 239 - 1; i >= 0; --i)
    {
        const uint8_t feedback = (info239[i] & 1u) ^ reg[15];
        for (int j = 15; j > 0; --j)
            reg[j] = static_cast<uint8_t>( reg[j - 1] ^ (G_COEFFS[j] & feedback) );
        reg[0] = static_cast<uint8_t>( G_COEFFS[0] & feedback ); // g0=1 => reg[0]=feedback
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

std::array<uint8_t,256> bch_255_239_encode(const std::vector<uint8_t>& info239)
{
    std::array<uint8_t,256> out{};
    // 拷入 239 个信息位（不足补 0）
    for (int i = 0; i < 239; ++i)
        out[i] = (i < static_cast<int>(info239.size())) ? (info239[i] & 1u) : 0u;

    // 计算 16 位校验并写到 [239..254]
    const auto par = bch_255_239_parity(info239);
    for (int j = 0; j < 16; ++j)
        out[239 + j] = par[j];

    // 计算整体偶校验位（让 256 位异或和为 0）
    uint8_t acc = 0;
    for (int i = 0; i < 255; ++i) acc ^= out[i];
    out[255] = acc;

    return out;
}

} // namespace newcode
