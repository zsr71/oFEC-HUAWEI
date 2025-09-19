#include "newcode/qam_llr.hpp"
#include "newcode/qam.hpp"   // 复用你的映射，保证一致性
#include "newcode/awgn.hpp"  // 用 ebn0_to_sigma
#include <cmath>
#include <stdexcept>

namespace newcode {

// 用现有 qam_modulate 构造与解调一致的星座查找表 S[j]
static std::vector<std::complex<float>>
build_constellation(unsigned n_bps)
{
    if (n_bps == 0 || (n_bps & 1u))
        throw std::invalid_argument("qam demod: n_bps must be a positive even number.");

    const unsigned M = 1u << n_bps;
    std::vector<uint8_t> bits;
    bits.reserve(M * n_bps);
    for (unsigned j = 0; j < M; ++j)
        for (unsigned b = 0; b < n_bps; ++b)
            bits.push_back( (j >> b) & 1u ); // 与调制/解调的位序完全一致 (LSB=bit0)

    // qam_modulate 会按顺序每 n_bps 个比特映射 1 个符号
    auto S = qam_modulate(bits, n_bps);
    return S; // size=M
}

std::vector<float>
qam_llr_maxlog(const std::vector<std::complex<float>>& y,
               unsigned n_bps,
               float sigma)
{
    if (y.empty()) return {};
    if (n_bps == 0 || (n_bps & 1u))
        throw std::invalid_argument("qam_llr_maxlog: n_bps must be a positive even number.");
    if (!(sigma > 0.f))
        throw std::invalid_argument("qam_llr_maxlog: sigma must be > 0.");

    const float inv_sigma2 = 1.f / (2.f * sigma * sigma);
    const unsigned M = 1u << n_bps;

    static thread_local unsigned cached_n_bps = 0;
    static thread_local std::vector<std::complex<float>> cached_S;
    if (cached_n_bps != n_bps || cached_S.size() != M) {
        cached_S = build_constellation(n_bps);
        cached_n_bps = n_bps;
    }
    const auto& S = cached_S;

    const size_t N_bits = y.size() * n_bps;
    std::vector<float> LLR(N_bits);

    for (size_t n = 0; n < N_bits; ++n)
    {
        const unsigned b = (unsigned)(n % n_bps);     // 比特位置
        const size_t   k = n / n_bps;                 // 符号索引
        const auto     yk = y[k];

        float L0 = -1e30f, L1 = -1e30f;               // -inf 的安全近似
        for (unsigned j = 0; j < M; ++j)
        {
            const float d2 = std::norm(yk - S[j]);    // |y - s_j|^2
            const float met = - d2 * inv_sigma2;      // -||y-s||^2 / (2σ^2)
            if (((j >> b) & 1u) == 0u)  { if (met > L0) L0 = met; }
            else                        { if (met > L1) L1 = met; }
        }
        LLR[n] = L0 - L1; // >0 表示更偏向 bit=0
    }
    return LLR;
}

std::vector<float>
qam_llr_from_ebn0(const std::vector<std::complex<float>>& y,
                  unsigned n_bps,
                  float ebn0_dB,
                  float code_rate)
{
    const float sigma = ebn0_to_sigma(ebn0_dB, n_bps, code_rate);
    return qam_llr_maxlog(y, n_bps, sigma);
}

} // namespace newcode
