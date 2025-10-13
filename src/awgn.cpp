#include "newcode/awgn.hpp"
#include <random>
#include <cmath>
#include <stdexcept>

namespace newcode {

float ebn0_to_sigma(float ebn0_dB, unsigned bits_per_symbol, float code_rate)
{
    if (bits_per_symbol == 0) throw std::invalid_argument("bits_per_symbol must be > 0");
    if (code_rate <= 0.f || code_rate > 1.f) throw std::invalid_argument("code_rate must be in (0,1]");

    const float ebn0_lin = std::pow(10.0f, ebn0_dB / 10.0f);
    const float esn0_lin = ebn0_lin * (bits_per_symbol);
    const float N0       = 1.0f / esn0_lin;      // 因为 Es=1
    const float sigma    = std::sqrt(N0 * 0.5f); // 每一维方差 N0/2
    return sigma;
}

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
