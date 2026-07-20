// path: newcode/qam_llr.cpp
#include "newcode/rx/demod/qam_llr.hpp"
#include "newcode/tx/mod/qam.hpp"
#include "newcode/channel/awgn.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include <complex>

namespace demod {

// 与原来一致：用现有映射构造星座
static std::vector<std::complex<float>>
build_constellation(unsigned n_bps)
{
    if (n_bps == 0)
        throw std::invalid_argument("qam demod: n_bps must be > 0.");
    if (n_bps > 1 && (n_bps & 1u))
        throw std::invalid_argument("qam demod: n_bps must be 1 or a positive even number.");

    const unsigned M = 1u << n_bps;
    std::vector<uint8_t> bits;
    bits.reserve(M * n_bps);
    for (unsigned j = 0; j < M; ++j)
        for (unsigned b = 0; b < n_bps; ++b)
            bits.push_back( (j >> b) & 1u ); // LSB=bit0

    return mod::qam_modulate(bits, n_bps); // size=M
}

// ===== 精确（log-sum-exp）版 =====
std::vector<float>
qam_llr_logsumexp(const std::vector<std::complex<float>>& y,
                  unsigned n_bps,
                  float sigma)
{
    if (y.empty()) return {};
    if (!(sigma > 0.f))
        throw std::invalid_argument("qam_llr_logsumexp: sigma must be > 0.");

    if (n_bps == 0)
        throw std::invalid_argument("qam_llr_logsumexp: n_bps must be > 0.");

    if (n_bps == 1)
    {
        const float inv_sigma_sq = 1.0f / (sigma * sigma);
        std::vector<float> LLR(y.size());
        for (size_t i = 0; i < y.size(); ++i)
        {
            //const float llr = 2.0f * y[i].real() * inv_sigma_sq;
            const float llr = y[i].real();
            LLR[i] = llr;
        }
        return LLR;
    }

    if (n_bps & 1u)
        throw std::invalid_argument("qam_llr_logsumexp: n_bps must be 1 or a positive even number.");

    const float inv_sigma2 = 1.f / (2.f * sigma * sigma); // 1/N0
    const unsigned M = 1u << n_bps;

    // 线程本地缓存星座
    static thread_local unsigned cached_n_bps = 0;
    static thread_local std::vector<std::complex<float>> cached_S;
    if (cached_n_bps != n_bps || cached_S.size() != M) {
        cached_S = build_constellation(n_bps);
        cached_n_bps = n_bps;
    }
    const auto& S = cached_S;

    const size_t N_bits = y.size() * n_bps;
    std::vector<float> LLR(N_bits);

    const float neg_inf = -std::numeric_limits<float>::infinity();

    for (size_t n = 0; n < N_bits; ++n)
    {
        const unsigned b = static_cast<unsigned>(n % n_bps); // 比特位索引（符号内）
        const size_t   k = n / n_bps;                        // 符号索引
        const auto     yk = y[k];

        // 第1遍：分别找 bit=0 / bit=1 的最大度量（避免溢出）
        float m0 = neg_inf, m1 = neg_inf;
        for (unsigned j = 0; j < M; ++j) {
            const float d2  = std::norm(yk - S[j]);          // |y - s_j|^2
            const float met = - d2 * inv_sigma2;             // -||y-s||^2 / N0
            if (((j >> b) & 1u) == 0u) { if (met > m0) m0 = met; }
            else                        { if (met > m1) m1 = met; }
        }

        // 第2遍：做 log-sum-exp
        double sum0 = 0.0, sum1 = 0.0; // 用 double 累加更稳
        for (unsigned j = 0; j < M; ++j) {
            const float d2  = std::norm(yk - S[j]);
            const float met = - d2 * inv_sigma2;
            if (((j >> b) & 1u) == 0u) sum0 += std::exp(static_cast<double>(met - m0));
            else                        sum1 += std::exp(static_cast<double>(met - m1));
        }

        // LLR = logsumexp0 - logsumexp1
        const double lse0 = static_cast<double>(m0) + std::log(sum0);
        const double lse1 = static_cast<double>(m1) + std::log(sum1);
        LLR[n] = static_cast<float>((lse0 - lse1)) ;
    }
    return LLR;
}

// 便捷入口：从 Eb/N0 计算 sigma，再走精确版
std::vector<float>
qam_llr_from_ebn0(const std::vector<std::complex<float>>& y,
                  unsigned n_bps,
                  float ebn0_dB,
                  float code_rate)
{
    const float sigma = channel::ebn0_to_sigma(ebn0_dB, n_bps, code_rate);
    return qam_llr_logsumexp(y, n_bps, sigma);
}

} // namespace newcode
