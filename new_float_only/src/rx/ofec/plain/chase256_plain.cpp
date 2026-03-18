#include "new_float_only/params.hpp"
#include "new_float_only/rx/ofec/chase/chase256.hpp"
#include "new_float_only/common/bch/bch_255_239.hpp"

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
#include "detail/chase256_plain_params.ipp"
#include "detail/chase256_plain_reliability.ipp"
#include "detail/chase256_plain_patterns.ipp"
#include "detail/chase256_plain_parity.ipp"
#include "detail/chase256_plain_trace_csv.ipp"
#include "detail/chase256_plain_impl.ipp"

namespace chase {

} // namespace chase
