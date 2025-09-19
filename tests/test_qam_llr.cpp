#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include "newcode/qam_llr.hpp"
#include <random>
#include <cassert>
#include <cmath>
#include <iostream>

static inline uint8_t hard_bit(float llr) { return llr < 0.f ? 1u : 0u; }

int main()
{
    using namespace newcode;

    // 场景：QPSK（n_bps=2）、R=1，目标 Eb/N0 = 6 dB
    const unsigned n_bps   = 2;
    const float    code_R  = 1.0f;
    const float    EbN0dB  = 6.0f;
    const uint32_t seed    = 4242;

    // 生成比特并调制
    std::mt19937 rng(2025);
    std::bernoulli_distribution bd(0.5);
    const size_t n_syms = 20000;
    std::vector<uint8_t> bits(n_syms * n_bps);
    for (auto &b : bits) b = (uint8_t)bd(rng);

    auto x = qam_modulate(bits, n_bps);

    // 加噪
    auto y = add_awgn(x, EbN0dB, n_bps, code_R, seed);

    // 解调 LLR（用同样的 Eb/N0->sigma 换算）
    const float sigma = ebn0_to_sigma(EbN0dB, n_bps, code_R);
    auto llr = qam_llr_maxlog(y, n_bps, sigma);

    // 评估 BER（硬判决）
    size_t err = 0;
    for (size_t i = 0; i < bits.size(); ++i)
        if (hard_bit(llr[i]) != (bits[i] & 1u)) err++;

    const double BER = (double)err / (double)bits.size();
    std::cout << "QAM LLR test: Eb/N0=" << EbN0dB << " dB, BER=" << BER << "\n";

    // 宽松阈值（6 dB 下 QPSK BER 应该远小于 0.1）
    assert(BER < 0.1);

    // 形状检查
    assert(llr.size() == bits.size());

    return 0;
}
