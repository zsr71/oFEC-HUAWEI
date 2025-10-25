#include "newcode/bch_255_239.hpp"
#include <cassert>
#include <random>
#include <iostream>
#include <array>
#include <vector>
#include <iomanip>
#include <string>
#include <cstring>
#include <algorithm>

// ----------------------------- 与实现一致的 G_COEFFS（用于 LFSR 伴随式） -----------------------------
static constexpr uint8_t G_COEFFS[16] = {
    1,1,0,0,0,1,1,0, 1,1,1,1,0,1,1,0
};

static constexpr int N_TOTAL = 256;
static constexpr int N_CORE  = 255;
static constexpr int PAR_IDX = 255;

// ----------------------------- LFSR 余数 / 伴随式（16 位） -----------------------------
static std::array<uint8_t,16>
lfsr_remainder_from_core(const std::array<uint8_t,255>& cw255)
{
    // 约定：cw255 = [info(239), parity(16)]，与编码保持系统码布局一致
    std::array<uint8_t,16> reg{}; // 全 0
    // 先喂 239 个信息位（MSB->LSB：i=238..0）
    for (int i = 239 - 1; i >= 0; --i)
    {
        const uint8_t feedback = (cw255[i] & 1u) ^ reg[15];
        for (int j = 15; j > 0; --j)
            reg[j] = static_cast<uint8_t>( reg[j - 1] ^ (G_COEFFS[j] & feedback) );
        reg[0] = static_cast<uint8_t>( G_COEFFS[0] & feedback );
    }
    // 再喂 16 位 parity（MSB->LSB：i=15..0）
    for (int i = 15; i >= 0; --i)
    {
        const uint8_t feedback = (cw255[239 + i] & 1u) ^ reg[15];
        for (int j = 15; j > 0; --j)
            reg[j] = static_cast<uint8_t>( reg[j - 1] ^ (G_COEFFS[j] & feedback) );
        reg[0] = static_cast<uint8_t>( G_COEFFS[0] & feedback );
    }
    return reg;
}

static uint16_t syndrome16_from_core(const std::array<uint8_t,255>& cw255)
{
    auto reg = lfsr_remainder_from_core(cw255);
    uint16_t syn = 0;
    for (int i = 15; i >= 0; --i) syn = (uint16_t)((syn << 1) | (reg[i] & 1u));
    return syn;
}

static bool overall_even_parity_256(const std::array<uint8_t,256>& cw)
{
    uint8_t acc = 0;
    for (int i = 0; i < 256; ++i) acc ^= (cw[i] & 1u);
    return acc == 0;
}

static uint8_t parity256_from255(const std::array<uint8_t,255>& cw255)
{
    uint8_t acc = 0;
    for (int i = 0; i < 255; ++i) acc ^= cw255[i];
    return acc;
}

// ----------------------------- GF(256) S1..S4（与译码器相同的表法） -----------------------------
namespace gfmini {
struct GF256
{
    static constexpr int PRIM = 0x11D; // x^8 + x^4 + x^3 + x^2 + 1
    static constexpr int N    = 255;

    static bool      inited;
    static uint8_t   alpha_to[N];
    static int16_t   index_of[256];

    static void init()
    {
        if (inited) return;
        index_of[0] = -1;
        uint16_t a = 1;
        for (int i = 0; i < N; ++i)
        {
            alpha_to[i] = static_cast<uint8_t>(a);
            index_of[alpha_to[i]] = static_cast<int16_t>(i);
            a <<= 1;
            if (a & 0x100) a ^= PRIM;
        }
        inited = true;
    }
};
bool GF256::inited = false;
uint8_t GF256::alpha_to[GF256::N];
int16_t GF256::index_of[256];

static void syndromes_S1_S4(const std::array<uint8_t,255>& cw255, uint8_t S[4])
{
    GF256::init();
    S[0]=S[1]=S[2]=S[3]=0;
    for (int j = 0; j < GF256::N; ++j)
    {
        if ((cw255[j] & 1u) == 0) continue;
        int e1 =  j                % GF256::N;
        int e2 = (j + e1)          % GF256::N;   // 2*j
        int e3 = (j + e2)          % GF256::N;   // 3*j
        int e4 = (j + e3)          % GF256::N;   // 4*j
        S[0] ^= GF256::alpha_to[e1];
        S[1] ^= GF256::alpha_to[e2];
        S[2] ^= GF256::alpha_to[e3];
        S[3] ^= GF256::alpha_to[e4];
    }
}
} // namespace gfmini

// ----------------------------- 打印工具 -----------------------------
static std::string bits_to_string(const std::vector<uint8_t>& v, int max_head=32, int max_tail=32)
{
    std::string s;
    const int n = (int)v.size();
    if (n <= max_head + max_tail)
    {
        for (int i = 0; i < n; ++i) s += char('0' + (v[i] & 1u));
    }
    else
    {
        for (int i = 0; i < max_head; ++i) s += char('0' + (v[i] & 1u));
        s += "...";
        for (int i = n - max_tail; i < n; ++i) s += char('0' + (v[i] & 1u));
        s += "  (len=" + std::to_string(n) + ")";
    }
    return s;
}

template<size_t N>
static std::string bits_to_string(const std::array<uint8_t,N>& a, int max_head=48, int max_tail=48)
{
    std::vector<uint8_t> v(a.begin(), a.end());
    return bits_to_string(v, max_head, max_tail);
}

static void print_syndrome16(uint16_t syn)
{
    std::cout << "LFSR Syndrome16 = 0x" << std::hex << std::setw(4) << std::setfill('0')
              << syn << std::dec << "  (";
    for (int b = 15; b >= 0; --b) std::cout << ((syn >> b) & 1);
    std::cout << ")\n";
}

static void print_S1S4(const uint8_t S[4])
{
    std::cout << "GF(256) S1..S4 = {";
    for (int i=0;i<4;++i)
    {
        if (i) std::cout << ", ";
        std::cout << "0x" << std::hex << std::setw(2) << std::setfill('0')
                  << (int)S[i] << std::dec;
    }
    std::cout << "}\n";
}

static int hamming_diff_255(const std::array<uint8_t,255>& a, const std::array<uint8_t,255>& b)
{
    int e=0; for (int i=0;i<255;++i) e += (a[i]^b[i]) & 1u; return e;
}

// ----------------------------- 一个场景：注错 -> 调你的译码器 -> 打印 -----------------------------
static void run_case(const char* title,
                     const std::array<uint8_t,255>& cw255_true,
                     const std::vector<int>& flip_pos) // 仅 0..254
{
    using namespace newcode;

    std::array<uint8_t,255> noisy = cw255_true;
    for (int idx : flip_pos)
        if (0 <= idx && idx < 255) noisy[idx] ^= 1u;

    // 伴随式 / 症状（噪声前）
    uint16_t syn16_noisy = syndrome16_from_core(noisy);
    uint8_t  S_noisy[4];  gfmini::syndromes_S1_S4(noisy, S_noisy);

    // 译码
    std::array<uint8_t,255> decoded{};
    bool ok = bch_255_239_decode_hiho_cw_255(noisy.data(), decoded.data());

    // 伴随式 / 症状（译码后）
    uint16_t syn16_dec = syndrome16_from_core(decoded);
    uint8_t  S_dec[4];   gfmini::syndromes_S1_S4(decoded, S_dec);

    // 结果对比
    int pre_err  = hamming_diff_255(noisy,   cw255_true);
    int post_err = hamming_diff_255(decoded, cw255_true);

    // 256 位偶校验（把 255 位扩展）
    std::array<uint8_t,256> noisy256{}, dec256{}, true256{};
    for (int i=0;i<255;++i) { noisy256[i]=noisy[i]; dec256[i]=decoded[i]; true256[i]=cw255_true[i]; }
    noisy256[255] = parity256_from255(noisy);
    dec256  [255] = parity256_from255(decoded);
    true256 [255] = parity256_from255(cw255_true);

    // 打印
    std::cout << "\n=== " << title << " ===\n";
    std::cout << "注入错误位置: ";
    if (flip_pos.empty()) std::cout << "(无)";
    for (size_t k=0;k<flip_pos.size();++k) std::cout << (k?", ":"") << flip_pos[k];
    std::cout << "\n";

    std::cout << "核心255位 Hamming:  注错前=" << pre_err
              << "  译码后=" << post_err
              << "  (译码返回=" << (ok ? "true" : "false") << ")\n";

    std::cout << "噪声码字：";
    print_syndrome16(syn16_noisy);
    print_S1S4(S_noisy);
    std::cout << "译码结果：";
    print_syndrome16(syn16_dec);
    print_S1S4(S_dec);

    std::cout << "256位整体偶校验： 噪声=" << (overall_even_parity_256(noisy256) ? "Even" : "Odd")
              << "  译码后="       << (overall_even_parity_256(dec256)   ? "Even" : "Odd")
              << "  真值="         << (overall_even_parity_256(true256)  ? "Even" : "Odd")
              << "\n";

    if (decoded == cw255_true && syn16_dec == 0 &&
        (S_dec[0]|S_dec[1]|S_dec[2]|S_dec[3]) == 0)
        std::cout << "⇒ 恢复为合法码字 ✅\n";
    else
        std::cout << "⇒ 未完全恢复为合法码字 ❌\n";
}

// ----------------------------- 主程序 -----------------------------
int main()
{
    using namespace newcode;

    // 1) 构造一帧可回归的数据并编码（系统码：[info239 | parity16]，再补整体偶校验位）
    std::mt19937 rng(12345);
    std::bernoulli_distribution bd(0.5);

    std::vector<uint8_t> info(239, 0);
    for (int i = 0; i < 239; ++i) info[i] = (uint8_t)bd(rng);

    auto par = bch_255_239_parity(info);  // 16 bits
    auto cw256 = bch_255_239_encode(info);// 256 bits

    // 拆出核心 255 位
    std::array<uint8_t,255> cw255{};
    for (int i=0;i<255;++i) cw255[i] = cw256[i];

    // 自检（应为合法码字）
    assert(syndrome16_from_core(cw255) == 0);
    {
        uint8_t S0[4]; gfmini::syndromes_S1_S4(cw255, S0);
        assert((S0[0]|S0[1]|S0[2]|S0[3]) == 0);
    }
    assert(overall_even_parity_256(cw256));

    // 打印编码前/后基础信息
    std::cout << "=== 原始信息(239) [head..tail] ===\n"
              << bits_to_string(info) << "\n\n";
    std::cout << "=== 16位 BCH 校验(239..254) ===\n"
              << bits_to_string(par, 16, 0) << "\n\n";
    std::cout << "=== 编码后码字(256) [head..tail] ===\n"
              << bits_to_string(cw256) << "\n";
    std::cout << "整体偶校验(256位) = " << (overall_even_parity_256(cw256) ? "Even/OK" : "Odd/FAIL") << "\n";
    std::cout << "真码字伴随式："; print_syndrome16(syndrome16_from_core(cw255));
    {
        uint8_t S0[4]; gfmini::syndromes_S1_S4(cw255, S0);
        std::cout << "真码字 "; print_S1S4(S0);
    }

    // 2) 场景：1/2/3 比特错（仅在 0..254 范围内注错）
    std::uniform_int_distribution<int> uni255(0, 254);
    int e1 = uni255(rng);
    int e2 = uni255(rng); while (e2 == e1) e2 = uni255(rng);
    int e3 = uni255(rng); while (e3 == e1 || e3 == e2) e3 = uni255(rng);

    run_case("场景A：1 比特错", cw255, {e1});
    run_case("场景B：2 比特错", cw255, {e1, e2});
    run_case("场景C：3 比特错（超过 t=2，预期大多无法完全恢复）", cw255, {e1, e2, e3});

    std::cout << "\nAll BCH(255,239) tests finished.\n";
    return 0;
}
