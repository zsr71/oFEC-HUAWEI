#pragma once

#include <cstddef>

namespace newcode {

struct Params;

/**
 * Chase-Pyndiah style component decoder operating on 255-bit BCH blocks
 * (no overall parity extension). Accepts soft input LLRs and produces
 * extrinsic information for the same 255 bit positions.
 *
 * @param Lin255  Input LLRs (length 255)
 * @param extr255 Output extrinsic LLRs (length 255)
 * @param p       Decoder control parameters (CHASE_L, CHASE_NTEST, beta, ALPHA)
 */
void chase_decode_255(const float* Lin255, float* extr255, const Params& p);

} // namespace newcode
