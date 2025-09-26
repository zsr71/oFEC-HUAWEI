#include "newcode/ber.hpp"
#include <algorithm>
#include <iostream>

namespace newcode {

BerStats compute_ber(const std::vector<uint8_t>& ref_bits,
                     const std::vector<uint8_t>& rx_bits)
{
    const std::size_t L = std::min(ref_bits.size(), rx_bits.size());

    std::size_t err = 0;
    for (std::size_t i = 0; i < L; ++i) {
        err += (ref_bits[i] ^ rx_bits[i]) & 1u;
    }

    BerStats s;
    s.errors = err;
    s.total  = L;
    s.ber    = (L == 0) ? 0.0 : static_cast<double>(err) / static_cast<double>(L);
    return s;
}

BerStats compute_and_print_ber(const std::vector<uint8_t>& ref_bits,
                               const std::vector<uint8_t>& rx_bits,
                               const char* label)
{
    BerStats s = compute_ber(ref_bits, rx_bits);
    std::cout << "[RESULT] " << (label ? label : "BER")
              << " BER=" << s.ber
              << "  (errs=" << s.errors << " / " << s.total << " compared)\n";
    return s;
}

} // namespace newcode
