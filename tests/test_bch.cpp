#include "newcode/bch_255_239.hpp"
#include <cassert>
#include <random>
#include <iostream>
#include <array>

// 复用同一组 G_COEFFS（测试内独立定义，避免和实现细节强耦合）
static constexpr uint8_t G_COEFFS[16] = {
    1,1,0,0,0,1,1,0, 1,1,1,1,0,1,1,0
};

// 用 LFSR 校验一个系统码 [info(239) | parity(16)] 是否为合法码字：
// 把 info MSB→LSB，再把 parity MSB→LSB 依次送入寄存器，若最后 16 位全 0，说明合法。
static bool lfsr_remainder_zero(const std::vector<uint8_t>& info239,
                                const std::array<uint8_t,16>& parity16)
{
    std::array<uint8_t,16> reg{};
    // info：i=238..0
    for (int i = 239 - 1; i >= 0; --i)
    {
        const uint8_t inb = (i < (int)info239.size()) ? (info239[i] & 1u) : 0u;
        const uint8_t feedback = inb ^ reg[15];
        for (int j = 15; j > 0; --j)
            reg[j] = static_cast<uint8_t>( reg[j - 1] ^ (G_COEFFS[j] & feedback) );
        reg[0] = static_cast<uint8_t>( G_COEFFS[0] & feedback );
    }
    // parity：MSB→LSB（15..0）
    for (int i = 15; i >= 0; --i)
    {
        const uint8_t feedback = (parity16[i] & 1u) ^ reg[15];
        for (int j = 15; j > 0; --j)
            reg[j] = static_cast<uint8_t>( reg[j - 1] ^ (G_COEFFS[j] & feedback) );
        reg[0] = static_cast<uint8_t>( G_COEFFS[0] & feedback );
    }
    // 余数应为 0
    for (auto v : reg) if (v) return false;
    return true;
}

// 校验 256 位整体偶校验位：异或和为 0
static bool overall_even_parity(const std::array<uint8_t,256>& cw)
{
    uint8_t acc = 0;
    for (int i = 0; i < 256; ++i) acc ^= (cw[i] & 1u);
    return acc == 0;
}

int main()
{
    using namespace newcode;

    // --- 用例 1：全 0 向量 ---
    {
        std::vector<uint8_t> info(239, 0);
        auto par = bch_255_239_parity(info);
        for (auto v : par) assert(v == 0);

        auto cw = bch_255_239_encode(info);
        // 前 239 为 0
        for (int i = 0; i < 239; ++i) assert(cw[i] == 0);
        // 16 位校验为 0
        for (int i = 0; i < 16; ++i) assert(cw[239 + i] == 0);
        // 整体偶校验位为 0
        assert(cw[255] == 0);

        // LFSR 余数为 0
        assert(lfsr_remainder_zero(info, par));
        // 256 位整体偶校验
        assert(overall_even_parity(cw));
    }

    // --- 用例 2：确定性随机（可回归） ---
    {
        std::mt19937 rng(12345);
        std::bernoulli_distribution bd(0.5);
        std::vector<uint8_t> info(239, 0);
        for (int i = 0; i < 239; ++i) info[i] = (uint8_t)bd(rng);

        auto par = bch_255_239_parity(info);
        auto cw  = bch_255_239_encode(info);

        // 编码的 239 信息应原样保留
        for (int i = 0; i < 239; ++i) assert(cw[i] == info[i]);
        // 编码尾部 16 位应等于 parity
        for (int i = 0; i < 16;  ++i) assert(cw[239 + i] == par[i]);

        // LFSR 余数为 0（系统码合法）
        assert(lfsr_remainder_zero(info, par));
        // 256 位整体偶校验
        assert(overall_even_parity(cw));
    }

    std::cout << "All BCH(255,239) tests passed.\n";
    return 0;
}
