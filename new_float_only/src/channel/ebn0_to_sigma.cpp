#include "new_float_only/channel/awgn.hpp"
#include <cmath>
#include <stdexcept>

namespace channel {

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

} // namespace new_float_only
