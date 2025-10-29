#pragma once
#include <vector>
#include <complex>
#include <cstdint>

namespace newcode {

// Eb/N0(dB) → 每一维高斯噪声标准差 sigma（假设调制已归一化到 Es=1）
float ebn0_to_sigma(float ebn0_dB, unsigned bits_per_symbol, float code_rate = 1.0f);

// 给复符号序列加 AWGN（按 Eb/N0 设置噪声）
// 必须从外部传入随机数种子 seed（相同 seed 复现实验）
std::vector<std::complex<float>>
add_awgn(const std::vector<std::complex<float>>& x,
         float ebn0_dB,
         unsigned bits_per_symbol,
         uint32_t seed);

} // namespace newcode
