// chase256_plain.cpp — Pyndiah '98 SISO (eqs. (14)–(17), (20), (21))
// -------------------------------------------------------------------------------------
// This file implements the soft-input/soft-output (SISO) Chase component decoder
// exactly following Pyndiah’s derivation:
//   • (14)–(17): soft output Λ_j = y_j + ω_j, where ω_j is computed from the
//     “best ML codeword” and the “best competing codeword with flipped bit j”
//     using correlation (inner-product) metrics equivalent to Euclidean criteria.
//   • (20): when no competing codeword exists for bit j, use a constant-reliability
//     fallback L0 whose magnitude reflects the average reliability; sign is that of
//     the ML decision at bit j.
//   • (21): the decoder must output ONLY EXTRINSIC information ω_j; the caller
//     shall form the next input as y(next) = y(channel) + α·ω, with α being a schedule.
// -------------------------------------------------------------------------------------
#include "newcode/chase256.hpp"
#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <cctype>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <type_traits>
#include <array>

// Plain-specific helper splits (kept as .ipp to avoid template link issues).
#include "detail/chase256_plain_constants.ipp"
#include "detail/chase256_plain_llr_adapters.ipp"
#include "detail/chase256_plain_params.ipp"
#include "detail/chase256_plain_reliability.ipp"
#include "detail/chase256_plain_patterns.ipp"
#include "detail/chase256_plain_parity.ipp"
#include "detail/chase256_plain_trace_csv.ipp"
#include "detail/chase256_plain_impl.ipp"

namespace newcode {

// ======================== explicit instantiations ========================
template void chase_decode_256_plain<float >(const float*,  const float*,  float*,  const Params&);
template void chase_decode_256_plain<int8_t>(const int8_t*, const int8_t*, float*, const Params&);
template void chase_decode_256_plain<float >(const float*,  float*,  const Params&);
template void chase_decode_256_plain<int8_t>(const int8_t*, float*, const Params&);

#define INSTANTIATE_CHASE256_PLAIN_QFLOAT(N) \
template void chase_decode_256_plain<newcode::qfloat<N>>( \
    const newcode::qfloat<N>*, const newcode::qfloat<N>*, float*, const Params&); \
template void chase_decode_256_plain<newcode::qfloat<N>>( \
    const newcode::qfloat<N>*, float*, const Params&);

INSTANTIATE_CHASE256_PLAIN_QFLOAT(2)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(3)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(4)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(5)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(6)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(7)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(8)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(9)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(10)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(11)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(12)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(13)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(14)
INSTANTIATE_CHASE256_PLAIN_QFLOAT(15)

#undef INSTANTIATE_CHASE256_PLAIN_QFLOAT

} // namespace newcode

