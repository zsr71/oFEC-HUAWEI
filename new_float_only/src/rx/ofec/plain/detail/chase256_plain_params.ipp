#pragma once

namespace new_float_only {
namespace detail {

// ----- select coefficients from new_float_only::Params -----
static inline void pick_cp(const new_float_only::Params& p, float& beta, float& alpha)
{
    // p.beta is reused as the fallback magnitude scale for L0 in (20)
    beta = p.beta;
    alpha = p.ALPHA;
}

} // namespace detail
} // namespace new_float_only

