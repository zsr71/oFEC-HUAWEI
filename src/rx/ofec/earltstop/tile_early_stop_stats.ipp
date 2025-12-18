#pragma once

namespace newcode {

template <typename LLR>
TileEarlyStopResult tile_early_stop_stats(const Matrix<LLR>& lin_matrix)
{
  TileEarlyStopResult res{};

  // Expected columns: 256 = 128 + 111 + 16(BCH parity) + 1(overall parity).
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  res.rows_total = rows;
  res.row_passed_flags.assign(rows, false);
  if (cols < 256) return res; // Defensive: unexpected shape, do not early-stop.
  std::array<uint8_t, 255> hard255{};

  bool all_pass = true;
  for (size_t r = 0; r < rows; ++r) {
    bool row_passed = true;
    for (int j = 0; j < 255; ++j) {
      const float v = llr_to_float(lin_matrix[r][static_cast<size_t>(j)]);
      hard255[static_cast<size_t>(j)] = (v < 0.0f) ? 1u : 0u;
    }

    if (!bch_255_239_syndromes_zero_cw_255(hard255.data())) {
      row_passed = false;
    } else {
      // 条件B：整体奇偶一致（extended parity）
      uint8_t parity255 = 0u;
      for (int j = 0; j < 255; ++j)
        parity255 ^= (hard255[static_cast<size_t>(j)] & 1u);

      const uint8_t overall =
          (llr_to_float(lin_matrix[r][255]) < 0.0f) ? 1u : 0u;

      if ((parity255 ^ overall) != 0u) {
        row_passed = false;
      }
    }

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
