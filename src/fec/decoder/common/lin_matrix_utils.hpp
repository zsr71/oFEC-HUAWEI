#pragma once

#include "newcode/bch_255_239.hpp"
#include "newcode/llr_utils.hpp"
#include "newcode/matrix.hpp"
#include "newcode/qfloat.hpp"

#include <array>
#include <cmath>
#include <cstdint>

namespace newcode {

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix)
{
  // 预期列数为 256 = 128(旧) + 111(新) + 16(BCH校验) + 1(整体奇偶)
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  if (cols < 256) return false; // 保守：尺寸异常则不早停

  std::array<uint8_t, 255> hard255{};
  std::array<uint8_t, 255> decoded255{};

  for (size_t r = 0; r < rows; ++r) {
    for (int j = 0; j < 255; ++j) {
      const float v = llr_to_float(lin_matrix[r][static_cast<size_t>(j)]);
      hard255[static_cast<size_t>(j)] = (v < 0.0f) ? 1u : 0u;
    }

    if (!bch_255_239_decode_hiho_cw_255(hard255.data(), decoded255.data()))
      return false;

    uint8_t parity255 = 0u;
    for (int j = 0; j < 255; ++j)
      parity255 ^= (hard255[static_cast<size_t>(j)] & 1u);
    const uint8_t overall = (llr_to_float(lin_matrix[r][255]) < 0.0f) ? 1u : 0u;

    if ((parity255 ^ overall) != 0u)
      return false;
  }

  return true;
}

template <typename LLR, typename Enable = void>
struct LinMatrixAdapter {
  using core_type = LLR;

  static core_type combine(const LLR& Lch, const LLR& La)
  {
    const float sum = llr_to_float(Lch) + llr_to_float(La);
    return llr_from_float<core_type>(sum);
  }

  static core_type channel(const LLR& v)
  {
    return llr_from_float<core_type>(llr_to_float(v));
  }
};

template <int NBITS, typename Store>
struct LinMatrixAdapter<qfloat<NBITS, Store>> {
  using core_type = float;

  static core_type combine(const qfloat<NBITS, Store>& Lch,
                           const qfloat<NBITS, Store>& La)
  {
    return static_cast<float>(Lch.code() + La.code());
  }

  static core_type channel(const qfloat<NBITS, Store>& v)
  {
    return static_cast<float>(v.code());
  }
};

template <typename LLR, typename Enable = void>
struct ExtrinsicQuantizer {
  static float quantize(float value) { return value; }
};

template <int NBITS, typename Store>
struct ExtrinsicQuantizer<qfloat<NBITS, Store>> {
  static float quantize(float value)
  {
    int code = static_cast<int>(std::lrint(value));
    const int lo = qfloat<NBITS, Store>::LO();
    const int hi = qfloat<NBITS, Store>::HI();
    if (code < lo) code = lo;
    if (code > hi) code = hi;
    qfloat<NBITS, Store> q;
    q.set_code(code);
    return llr_to_float(q);
  }
};

} // namespace newcode
