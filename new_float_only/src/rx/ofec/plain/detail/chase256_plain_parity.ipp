#pragma once

namespace new_float_only {
namespace detail {

// ----- overall (even) parity from 255-bit core (extend to 256 bits) -----
inline uint8_t parity256_from255(const uint8_t* cw255) {
    uint8_t acc = 0;
    for (int i = 0; i < BCH_N_CORE; ++i) acc ^= cw255[i];
    return acc;
}

} // namespace detail
} // namespace new_float_only

