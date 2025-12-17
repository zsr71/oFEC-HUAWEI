#pragma once

#include "newcode/bch_255_239.hpp"
#include "newcode/llr_utils.hpp"
#include "newcode/matrix.hpp"
#include "newcode/qfloat.hpp"

#include <array>
#include <cmath>
#include <cstdint>

namespace newcode {

struct TileEarlyStopResult {
  bool all_rows_passed = false;
  std::size_t rows_passed = 0;
  std::size_t rows_total = 0;
};

template <typename LLR>
TileEarlyStopResult tile_early_stop_stats(const Matrix<LLR>& lin_matrix)
{
  TileEarlyStopResult res{};

  // 棰勬湡鍒楁暟涓?256 = 128(鏃? + 111(鏂? + 16(BCH鏍￠獙) + 1(鏁翠綋濂囧伓)
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  res.rows_total = rows;
  if (cols < 256) return res; // 淇濆畧锛氬昂瀵稿紓甯稿垯涓嶆棭鍋?
  std::array<uint8_t, 255> hard255{};
  std::array<uint8_t, 255> decoded255{};

  bool all_pass = true;
  for (size_t r = 0; r < rows; ++r) {
    for (int j = 0; j < 255; ++j) {
      const float v = llr_to_float(lin_matrix[r][static_cast<size_t>(j)]);
      hard255[static_cast<size_t>(j)] = (v < 0.0f) ? 1u : 0u;
    }

    if (!bch_255_239_decode_hiho_cw_255(hard255.data(), decoded255.data())) {
      all_pass = false;
      continue;
    }

    uint8_t parity255 = 0u;
    for (int j = 0; j < 255; ++j)
      parity255 ^= (hard255[static_cast<size_t>(j)] & 1u);
    const uint8_t overall = (llr_to_float(lin_matrix[r][255]) < 0.0f) ? 1u : 0u;

    if ((parity255 ^ overall) != 0u) {
      all_pass = false;
      continue;
    }

    ++res.rows_passed;
  }

  res.all_rows_passed = all_pass && (rows > 0 ? res.rows_passed == rows : false);
  return res;
}

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix)
{
  return tile_early_stop_stats(lin_matrix).all_rows_passed;
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
