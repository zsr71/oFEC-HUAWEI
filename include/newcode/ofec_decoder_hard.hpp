#pragma once

#include <array>
#include <cstdint>

#include "newcode/params.hpp"
#include "newcode/common/qfloat/qfloat.hpp"

namespace newcode {

template <typename LLR>
bool perform_hard_decode(const std::array<LLR, 256>& Lin256,
                         const std::array<LLR, 256>& Lch256,
                         std::array<float, 256>& Y2,
                         const Params& p);

extern template bool perform_hard_decode<float>(const std::array<float, 256>&,
                                                const std::array<float, 256>&,
                                                std::array<float, 256>&,
                                                const Params&);
extern template bool perform_hard_decode<int8_t>(const std::array<int8_t, 256>&,
                                                 const std::array<int8_t, 256>&,
                                                 std::array<float, 256>&,
                                                 const Params&);

#define DECLARE_HARD_DECODE_QFLOAT(N) \
extern template bool perform_hard_decode<qfloat::qfloat<N>>( \
    const std::array<qfloat::qfloat<N>, 256>&, \
    const std::array<qfloat::qfloat<N>, 256>&, \
    std::array<float, 256>&, \
    const Params&);

DECLARE_HARD_DECODE_QFLOAT(2)
DECLARE_HARD_DECODE_QFLOAT(3)
DECLARE_HARD_DECODE_QFLOAT(4)
DECLARE_HARD_DECODE_QFLOAT(5)
DECLARE_HARD_DECODE_QFLOAT(6)
DECLARE_HARD_DECODE_QFLOAT(7)
DECLARE_HARD_DECODE_QFLOAT(8)
DECLARE_HARD_DECODE_QFLOAT(9)
DECLARE_HARD_DECODE_QFLOAT(10)
DECLARE_HARD_DECODE_QFLOAT(11)
DECLARE_HARD_DECODE_QFLOAT(12)
DECLARE_HARD_DECODE_QFLOAT(13)
DECLARE_HARD_DECODE_QFLOAT(14)
DECLARE_HARD_DECODE_QFLOAT(15)

#undef DECLARE_HARD_DECODE_QFLOAT

} // namespace newcode
