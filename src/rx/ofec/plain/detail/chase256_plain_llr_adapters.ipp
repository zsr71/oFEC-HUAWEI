#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

namespace newcode {
namespace detail {

// ----- LLR adapters -----
template<typename LLR>
inline float llr_to_float(LLR x) { return static_cast<float>(x); }

template<typename LLR, typename std::enable_if<std::is_floating_point<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return static_cast<LLR>(x); }

template<typename LLR, typename std::enable_if<std::is_integral<LLR>::value && std::is_signed<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x)
{
    const float lo = (float)std::numeric_limits<LLR>::min();
    const float hi = (float)std::numeric_limits<LLR>::max();
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return (LLR)std::lrintf(x);
}

template<typename LLR, typename std::enable_if<!std::is_arithmetic<LLR>::value, int>::type = 0>
inline LLR llr_from_float(float x) { return LLR::from_float(x); }

} // namespace detail
} // namespace newcode
