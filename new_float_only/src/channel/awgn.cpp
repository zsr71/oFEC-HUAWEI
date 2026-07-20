#include "new_float_only/channel/awgn.hpp"
#include <random>

namespace channel {

// 给复数基带符号序列添加 AWGN 噪声。
//
// 输入 x 被视为已经归一化好的发送符号序列。函数先根据 Eb/N0 和每个符号承载的
// 比特数计算噪声标准差 sigma，然后分别在 I/Q 两个维度上叠加独立高斯白噪声。
//
// 返回值 y 与输入 x 等长，表示通过 AWGN 信道后的接收复符号序列。
std::vector<std::complex<float>>
add_awgn(const std::vector<std::complex<float>>& x,
         float ebn0_dB,
         unsigned bits_per_symbol,
         uint32_t seed)
{
    // 把给定的 Eb/N0 换算成每个实维高斯噪声的标准差。
    const float sigma = ebn0_to_sigma(ebn0_dB, bits_per_symbol);

    // 随机数种子由外部传入，这样同一组参数下的信道实现可复现。
    std::mt19937 rng(seed);
    std::normal_distribution<float> gauss(0.0f, sigma);

    std::vector<std::complex<float>> y;
    y.reserve(x.size());
    for (auto s : x)
    {
        // 复数基带信号在 I/Q 两个维度上分别叠加独立同分布的高斯噪声。
        const float nI = gauss(rng);
        const float nQ = gauss(rng);
        y.emplace_back(s.real() + nI, s.imag() + nQ);
    }
    return y;
}

} // namespace new_float_only
