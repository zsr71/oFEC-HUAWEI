#include "newcode/bch_255_239.hpp"
#include <algorithm>
#include <vector>
#include <cmath>
#include <cstring>

namespace newcode
{
// ---------------- LFSR 奇偶（保留你的实现） ----------------
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

// ================== 硬判决译码（t=2，GF(2^8)） ==================
namespace {
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
    static inline uint8_t add(uint8_t x, uint8_t y) { return x ^ y; }
    static inline uint8_t mul(uint8_t x, uint8_t y)
    {
        if (!x || !y) return 0;
        int lx = index_of[x], ly = index_of[y];
        return alpha_to[(lx + ly) % N];
    }
    static inline uint8_t div(uint8_t x, uint8_t y)
    {
        if (!x) return 0;
        int lx = index_of[x], ly = index_of[y];
        return alpha_to[(lx - ly + N) % N];
    }
};

bool GF256::inited = false;
uint8_t GF256::alpha_to[GF256::N];
int16_t GF256::index_of[256];

static void compute_syndromes_1_4(const uint8_t* in255, uint8_t S[4])
{
    GF256::init();
    S[0]=S[1]=S[2]=S[3]=0;
    for (int j = 0; j < GF256::N; ++j)
    {
        if ((in255[j] & 1u) == 0) continue;
        int e1 =  j             % GF256::N;
        int e2 = (j + e1)       % GF256::N;   // 2*j
        int e3 = (j + e2)       % GF256::N;   // 3*j
        int e4 = (j + e3)       % GF256::N;   // 4*j
        S[0] ^= GF256::alpha_to[e1];
        S[1] ^= GF256::alpha_to[e2];
        S[2] ^= GF256::alpha_to[e3];
        S[3] ^= GF256::alpha_to[e4];
    }
}

// 轻量 BM（t=2）：σ(x)=1 + c1 x + c2 x^2
static int berlekamp_massey_t2(const uint8_t S[4], uint8_t sigma[3])
{
    using G = GF256;
    uint8_t C[3] = {1,0,0}, B[3] = {1,0,0};
    int L = 0, m = 1;
    uint8_t b = 1;

    for (int n = 0; n < 4; ++n)
    {
        uint8_t d = S[n];
        for (int i = 1; i <= L; ++i) d ^= G::mul(C[i], S[n-i]);

        if (d == 0) { m++; continue; }

        uint8_t T0=C[0], T1=C[1], T2=C[2];
        uint8_t db = G::div(d, b);
        if      (m == 1) { C[2] ^= G::mul(db, B[1]); C[1] ^= G::mul(db, B[0]); }
        else if (m >= 2) { C[2] ^= G::mul(db, B[0]); }

        if (2*L <= n) { L = n + 1 - L; B[0]=T0; B[1]=T1; B[2]=T2; b=d; m=1; }
        else          { m++; }
    }
    sigma[0]=C[0]; sigma[1]=C[1]; sigma[2]=C[2];
    return L; // 0/1/2
}

static int chien_and_correct(uint8_t* cw255, const uint8_t sigma[3], int L)
{
    using G = GF256;
    if (L == 0) return 0;

    int count = 0, locs[2] = {-1,-1};
    int reg1 = (sigma[1] ? G::index_of[sigma[1]] : -1);
    int reg2 = (L>=2 && sigma[2] ? G::index_of[sigma[2]] : -1);

    for (int i = 1; i <= G::N; ++i)
    {
        uint8_t q = 1;
        if (reg1 != -1) { reg1 = (reg1 + 1) % G::N; q ^= G::alpha_to[reg1]; }
        if (reg2 != -1) { reg2 = (reg2 + 2) % G::N; q ^= G::alpha_to[reg2]; }
        if (q == 0)
        {
            int pos = G::N - i; // 与 AFF3CT 一致：loc = 255 - i
            if (pos >= 0 && pos < G::N) { if (count < L) locs[count] = pos; count++; }
        }
    }
    if (count != L) return -1;
    for (int k = 0; k < L; ++k) if (locs[k] >= 0) cw255[locs[k]] ^= 1u;
    return L;
}
} // anon

bool bch_255_239_decode_hiho_cw_255(const uint8_t* in255, uint8_t* out255)
{
    using G = GF256;
    G::init();

    uint8_t cw[GF256::N];
    std::memcpy(cw, in255, GF256::N);

    uint8_t S[4]; compute_syndromes_1_4(cw, S);
    if ((S[0]|S[1]|S[2]|S[3]) == 0) { std::memcpy(out255, cw, GF256::N); return true; }

    uint8_t sigma[3]; int L = berlekamp_massey_t2(S, sigma);
    if (L < 0 || L > 2) { std::memcpy(out255, cw, GF256::N); return false; }

    int corr = chien_and_correct(cw, sigma, L);
    if (corr < 0) { std::memcpy(out255, in255, GF256::N); return false; }

    uint8_t S2[4]; compute_syndromes_1_4(cw, S2);
    bool ok = ((S2[0]|S2[1]|S2[2]|S2[3]) == 0);
    std::memcpy(out255, cw, GF256::N);
    return ok;
}

} // namespace newcode
