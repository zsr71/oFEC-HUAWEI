#pragma once
#include <array>
#include <vector>
#include <cstdint>
#include <cstddef>

namespace newcode
{
// ========== 生成多项式 g(y)（保持你原实现） ==========
static constexpr uint8_t G_COEFFS[16] = {
    /*y^0..y^15*/
    1, 1, 0, 0, 0, 1, 1, 0,
    1, 1, 1, 1, 0, 1, 1, 0
};

// ========== 编码（保持你原实现的 API） ==========
std::array<uint8_t,16>  bch_255_239_parity (const std::vector<uint8_t>& info239);
std::array<uint8_t,256> bch_255_239_encode (const std::vector<uint8_t>& info239);

// ========== 新增：硬判决译码（255 in / 255 out） ==========
// 适配 chase256：输入 255 位硬判决，输出 255 位纠正码字（不含整体奇偶）
// 返回 true=成功（合法或可纠错），false=失败
bool bch_255_239_decode_hiho_cw_255(const uint8_t* in255, uint8_t* out255);

// 便捷封装：若你传 256 位（丢弃第 256 位整体奇偶），对前 255 位译码
inline bool bch_255_239_decode_hiho_cw_256(const uint8_t* in256, uint8_t* out255)
{
    return bch_255_239_decode_hiho_cw_255(in256, out255);
}

} // namespace newcode
