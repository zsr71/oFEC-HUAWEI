#pragma once

#include <array>
#include <cstdint>

#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {

template <typename LLR>
bool perform_hard_decode(const std::array<LLR, 256>& Lin256,
                         const std::array<LLR, 256>& Lch256,
                         std::array<LLR, 256>& Y2,
                         const Params& p);

extern template bool perform_hard_decode<float>(const std::array<float, 256>&,
                                                const std::array<float, 256>&,
                                                std::array<float, 256>&,
                                                const Params&);
extern template bool perform_hard_decode<int8_t>(const std::array<int8_t, 256>&,
                                                 const std::array<int8_t, 256>&,
                                                 std::array<int8_t, 256>&,
                                                 const Params&);
extern template bool perform_hard_decode<qfloat<4>>(const std::array<qfloat<4>, 256>&,
                                                    const std::array<qfloat<4>, 256>&,
                                                    std::array<qfloat<4>, 256>&,
                                                    const Params&);
extern template bool perform_hard_decode<qfloat<5>>(const std::array<qfloat<5>, 256>&,
                                                    const std::array<qfloat<5>, 256>&,
                                                    std::array<qfloat<5>, 256>&,
                                                    const Params&);

} // namespace newcode
