#pragma once

#include "newcode/llr_utils.hpp"
#include "newcode/matrix.hpp"
#include "newcode/qfloat.hpp"

#include <cmath>

namespace newcode {

template <typename LLR, typename Enable = void>
struct LinMatrixAdapter {
  using core_type = LLR;

  static core_type combine(const LLR& Lch, const LLR& La);
  static core_type channel(const LLR& v);
};

template <int NBITS, typename Store>
struct LinMatrixAdapter<qfloat<NBITS, Store>> {
  using core_type = float;

  static core_type combine(const qfloat<NBITS, Store>& Lch,
                           const qfloat<NBITS, Store>& La);
  static core_type channel(const qfloat<NBITS, Store>& v);
};

template <typename LLR, typename Enable = void>
struct ExtrinsicQuantizer {
  static float quantize(float value);
};

template <int NBITS, typename Store>
struct ExtrinsicQuantizer<qfloat<NBITS, Store>> {
  static float quantize(float value);
};

} // namespace newcode

#include "ofec/common/lin_matrix_adapters.ipp"
