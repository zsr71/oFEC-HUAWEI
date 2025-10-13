#include "newcode/qam.hpp"
#include "newcode/awgn.hpp"
#include <cassert>
#include <random>
#include <iostream>
#include <cmath>

static inline bool approx(float a, float b, float tol) { return std::fabs(a - b) <= tol; }

int main()
{
    using namespace newcode;

    // 配置：QPSK（k=2），无编码 R=1
    const unsigned k = 2;
    const float    R = 1.0f;
    const float    EbN0_dB = 10.0f;
    const uint32_t seed1 = 12345;
    const uint32_t seed2 = 54321;

    // 生成比特并 QAM（你的 qam_modulate 会把 Es≈1）
    std::mt19937 rng_bits(2025);
    std::bernoulli_distribution bd(0.5);
    const size_t n_syms = 50000;
    std::vector<uint8_t> bits(n_syms * k);
    for (auto &b : bits) b = (uint8_t)bd(rng_bits);

    auto x = qam_modulate(bits, k);

    // 1) 可复现性：相同 seed → 相同噪声；不同 seed → 不同噪声
    auto y1 = add_awgn(x, EbN0_dB, k, seed1);
    auto y2 = add_awgn(x, EbN0_dB, k, seed1);
    auto y3 = add_awgn(x, EbN0_dB, k, seed2);

    assert(y1.size() == x.size() && y2.size() == x.size() && y3.size() == x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        assert(y1[i].real() == y2[i].real());
        assert(y1[i].imag() == y2[i].imag());
    }
    bool any_diff = false;
    for (size_t i = 0; i < x.size(); ++i) {
        if (y1[i] != y3[i]) { any_diff = true; break; }
    }
    assert(any_diff); // 不同 seed 应产生不同噪声轨迹

    // 2) 统计验证：估计 Eb/N0 与目标值接近；每维方差≈sigma^2
    auto y = y1; // 用 y1 做统计
    double sum_nI2 = 0.0, sum_nQ2 = 0.0, sum_n2 = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        const float nI = y[i].real() - x[i].real();
        const float nQ = y[i].imag() - x[i].imag();
        sum_nI2 += nI * nI;
        sum_nQ2 += nQ * nQ;
        sum_n2  += nI * nI + nQ * nQ;
    }
    const double var_I = sum_nI2 / x.size();
    const double var_Q = sum_nQ2 / x.size();
    const double N0_est = sum_n2 / x.size();       // 复噪声功率 ≈ N0
    const double EsN0_est_lin = 1.0 / N0_est;      // Es=1
    const double EbN0_est_lin = EsN0_est_lin / (k * R);
    const double EbN0_est_dB  = 10.0 * std::log10(EbN0_est_lin);

    const float sigma = ebn0_to_sigma(EbN0_dB, k, R);
    const float var_th = sigma * sigma;

    assert(approx((float)var_I, var_th, 0.005f));
    assert(approx((float)var_Q, var_th, 0.005f));
    assert(std::fabs(EbN0_est_dB - EbN0_dB) < 0.4);

    std::cout << "AWGN test passed. var_I=" << var_I
              << " var_Q=" << var_Q
              << " EbN0_est_dB=" << EbN0_est_dB << "\n";
    return 0;
}
