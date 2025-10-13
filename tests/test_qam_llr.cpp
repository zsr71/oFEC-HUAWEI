// path: tests/test_qam_llr.cpp
#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include "newcode/qam_llr.hpp"
#include <random>
#include <cassert>
#include <cmath>
#include <iostream>

static inline uint8_t hard_bit(float llr) { return llr < 0.f ? 1u : 0u; }

int main() {
    using namespace newcode;

    // 场景：QPSK（n_bps=2）、R=1，目标 Eb/N0 = 6 dB
    const unsigned n_bps   = 2;
    const float    code_R  = 1.0f;
    const float    EbN0dB  = 6.0f;
    const uint32_t seed    = 4242;

    // 生成随机比特并调制
    std::mt19937 rng(2025);
    std::bernoulli_distribution bd(0.5);
    const size_t n_syms = 20000;
    std::vector<uint8_t> bits(n_syms * n_bps);
    for (auto& b : bits) b = static_cast<uint8_t>(bd(rng));

    auto x = qam_modulate(bits, n_bps);

    // 加噪（注意：带上 code_rate）
    auto y = add_awgn(x, EbN0dB, n_bps, seed);

    // sigma ← Eb/N0 换算；N0=2*sigma^2
    const float sigma = ebn0_to_sigma(EbN0dB, n_bps, code_R);

    // 精确（log-sum-exp）LLR
    auto llr_exact = qam_llr_logsumexp(y, n_bps, sigma);



    // 评估 BER（硬判决）
    auto ber_from_llr = [&](const std::vector<float>& llr) {
        size_t err = 0;
        for (size_t i = 0; i < bits.size(); ++i)
            if (hard_bit(llr[i]) != (bits[i] & 1u)) err++;
        return static_cast<double>(err) / static_cast<double>(bits.size());
    };

    const double BER_exact  = ber_from_llr(llr_exact);

    std::cout << "QPSK Eb/N0=" << EbN0dB
              << " dB  BER_exact="  << BER_exact;

    // 理论上 QPSK 在 6 dB 处 BER ~ 2.4e-3，这里给个宽松但合理的阈值
    assert(BER_exact  < 0.01);

    // 形状检查
    assert(llr_exact.size()  == bits.size());

    return 0;
}
