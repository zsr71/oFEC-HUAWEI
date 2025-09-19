#include "newcode/qam.hpp"
#include <cassert>
#include <iostream>
#include <random>
#include <cmath>

static inline bool approx(float a, float b, float eps = 1e-6f)
{
    return std::fabs(a - b) <= eps;
}

int main()
{
    using newcode::qam_modulate;

    // --- 用例 1：4-QAM（QPSK）固定映射检查 ---
    {
        // 分组按 [b_I, b_Q]
        std::vector<uint8_t> bits = {0,0, 0,1, 1,0, 1,1};
        auto syms = qam_modulate(bits, 2);
        assert(syms.size() == 4);

        const float inv_sqrt2 = 1.0f / std::sqrt(2.0f);
        // 00 -> (+1,+1)/√2
        assert(approx(syms[0].real(), +inv_sqrt2));
        assert(approx(syms[0].imag(), +inv_sqrt2));
        // 01 -> (+1,-1)/√2
        assert(approx(syms[1].real(), +inv_sqrt2));
        assert(approx(syms[1].imag(), -inv_sqrt2));
        // 10 -> (-1,+1)/√2
        assert(approx(syms[2].real(), -inv_sqrt2));
        assert(approx(syms[2].imag(), +inv_sqrt2));
        // 11 -> (-1,-1)/√2
        assert(approx(syms[3].real(), -inv_sqrt2));
        assert(approx(syms[3].imag(), -inv_sqrt2));
    }

    // --- 用例 2：bits 不是 n_bps 的整数倍时补零，并检查长度 ---
    {
        std::vector<uint8_t> bits = {1, 0, 1}; // 3 位 → 2 个符号（补 0）
        auto syms = qam_modulate(bits, 2);
        assert(syms.size() == 2);
    }

    // --- 用例 3：16-QAM 统计能量接近 1（归一化） ---
    {
        std::mt19937 rng(1234);
        std::bernoulli_distribution bd(0.5);
        std::vector<uint8_t> bits(4 * 5000, 0); // 5000 个符号
        for (auto &b : bits) b = (uint8_t)bd(rng);

        auto syms = qam_modulate(bits, 4);
        double Es_avg = 0.0;
        for (auto s : syms) Es_avg += (double)(s.real()*s.real() + s.imag()*s.imag());
        Es_avg /= (double)syms.size();
        // 允许一点统计误差
        assert(std::fabs(Es_avg - 1.0) < 0.05);
    }

    std::cout << "All QAM tests passed.\n";
    return 0;
}
