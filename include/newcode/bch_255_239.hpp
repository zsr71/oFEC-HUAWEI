#pragma once
#include <array>
#include <vector>
#include <cstdint>

namespace newcode {

// 计算 16 位 BCH(255,239) 校验（信息位不足 239 会自动 0 填充）
std::array<uint8_t,16> bch_255_239_parity(const std::vector<uint8_t>& info239);

// 生成完整系统码（长度 256 = 239 info + 16 parity + 1 overall parity，顺序为 [info | parity | overall]）
std::array<uint8_t,256> bch_255_239_encode(const std::vector<uint8_t>& info239);

} // namespace newcode
