#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/common/qfloat/qfloat.hpp"
#include "newcode/params.hpp"
#include "newcode/rx/ofec/chase/chase256.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "detail/chase256_plain_constants.ipp"
#include "detail/chase256_plain_llr_adapters.ipp"
#include "detail/chase256_plain_params.ipp"
#include "detail/chase256_plain_reliability.ipp"
#include "detail/chase256_plain_patterns.ipp"
#include "detail/chase256_plain_parity.ipp"
#include "detail/chase256_plain_trace_csv.ipp"
#include "detail/chase256_global_pair_impl.ipp"

namespace chase {

template void chase_decode_256_global_pair<float >(const float*,  const float*,  float*,  const newcode::Params&);
template void chase_decode_256_global_pair<int8_t>(const int8_t*, const int8_t*, float*, const newcode::Params&);
template void chase_decode_256_global_pair<float >(const float*,  float*,  const newcode::Params&);
template void chase_decode_256_global_pair<int8_t>(const int8_t*, float*, const newcode::Params&);

#define INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(N) \
template void chase_decode_256_global_pair<qfloat::qfloat<N>>( \
    const qfloat::qfloat<N>*, const qfloat::qfloat<N>*, float*, const newcode::Params&); \
template void chase_decode_256_global_pair<qfloat::qfloat<N>>( \
    const qfloat::qfloat<N>*, float*, const newcode::Params&);

INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(2)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(3)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(4)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(5)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(6)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(7)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(8)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(9)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(10)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(11)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(12)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(13)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(14)
INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT(15)

#undef INSTANTIATE_CHASE256_GLOBAL_PAIR_QFLOAT

} // namespace chase
