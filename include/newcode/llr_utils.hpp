#pragma once

#include <cmath>
#include <limits>
#include <type_traits>

namespace newcode {

template <typename LLR>
inline float llr_to_float(LLR x)
{
  return static_cast<float>(x);
}

template <typename LLR>
inline LLR llr_from_float(float x)
{
  if constexpr (std::is_floating_point_v<LLR>) {
    return static_cast<LLR>(x);
  } else if constexpr (std::is_integral_v<LLR> && std::is_signed_v<LLR>) {
    const float lo = static_cast<float>(std::numeric_limits<LLR>::min());
    const float hi = static_cast<float>(std::numeric_limits<LLR>::max());
    if (x < lo) x = lo;
    if (x > hi) x = hi;
    return static_cast<LLR>(std::lrintf(x));
  } else {
    return LLR::from_float(x);
  }
}

} // namespace newcode

