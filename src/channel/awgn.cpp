#include "newcode/awgn.hpp"
#include <random>

namespace newcode {

std::vector<std::complex<float>>
add_awgn(const std::vector<std::complex<float>>& x,
         float ebn0_dB,
         unsigned bits_per_symbol,
         uint32_t seed)
{
    const float sigma = ebn0_to_sigma(ebn0_dB, bits_per_symbol);

    std::mt19937 rng(seed); // 必须由外部传入，确保可复现
    std::normal_distribution<float> gauss(0.0f, sigma);

    std::vector<std::complex<float>> y;
    y.reserve(x.size());
    for (auto s : x)
    {
        const float nI = gauss(rng);
        const float nQ = gauss(rng);
        y.emplace_back(s.real() + nI, s.imag() + nQ);
    }
    return y;
}

} // namespace newcode
