#pragma once

#include "newcode/common/bch/bch_255_239.hpp"

namespace newcode {

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop_v1(const matrix::Matrix<LLR>& lin_matrix,
                                              const Params& p)
{
  (void)p;
  TileEarlyStopResult res{};

  // Expected columns: 256 = 128 + 111 + 16(BCH parity) + 1(overall parity).
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  res.rows_total = rows;
  res.row_passed_flags.assign(rows, false);
  res.row_details.assign(rows, TileEarlyStopRowDetail{});
  if (cols < 256) return res; // Defensive: unexpected shape, do not early-stop.
  std::array<uint8_t, 255> hard255{};

  bool all_pass = true;
  for (size_t r = 0; r < rows; ++r) {
    auto& row_detail = res.row_details[r];
    for (int j = 0; j < 255; ++j) {
      const float v = qfloat::llr_to_float(lin_matrix[r][static_cast<size_t>(j)]);
      hard255[static_cast<size_t>(j)] = (v < 0.0f) ? 1u : 0u;
    }

    const uint8_t syndrome_bits =
        bch::bch_255_239_syndrome_nonzero_mask_cw_255(hard255.data());
    row_detail.syndrome_bits = syndrome_bits;
    const bool bch_passed = (syndrome_bits == 0u);

    uint8_t parity255 = 0u;
    for (int j = 0; j < 255; ++j)
      parity255 ^= (hard255[static_cast<size_t>(j)] & 1u);

    const uint8_t overall =
        (qfloat::llr_to_float(lin_matrix[r][255]) < 0.0f) ? 1u : 0u;
    const bool overall_passed = ((parity255 ^ overall) == 0u);

    const bool row_passed =
        (!p.EARLY_STOP_COND_V1_REQUIRE_BCH || bch_passed) &&
        (!p.EARLY_STOP_COND_V1_REQUIRE_OVERALL || overall_passed);
    row_detail.bch_passed = bch_passed;
    row_detail.overall_parity_passed = overall_passed;
    row_detail.row_passed = row_passed;

    if (row_passed) {
      ++res.rows_passed;
    } else {
      all_pass = false;
    }
    res.row_passed_flags[r] = row_passed;
  }

  res.all_rows_passed = all_pass && (rows > 0 ? res.rows_passed == rows : false);
  return res;
}

} // namespace newcode
