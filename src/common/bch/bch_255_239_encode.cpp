#include "newcode/bch_255_239.hpp"

namespace newcode
{

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
