#pragma once

namespace newcode {
namespace detail {

// ----- select coefficients from newcode::Params -----
static inline void pick_cp(const newcode::Params& p, float& beta, float& alpha)
{
    // p.beta is reused as the fallback magnitude scale for L0 in (20)
    beta = p.beta;
    alpha = p.ALPHA;
}

} // namespace detail
} // namespace newcode

