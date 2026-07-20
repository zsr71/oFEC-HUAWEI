#pragma once

#include <algorithm>

namespace newcode {

inline TileEarlyStopResult apply_group_bound_early_stop(
    const TileEarlyStopResult& raw_result,
    int bind_group_size) {
  if (bind_group_size <= 1 || raw_result.row_passed_flags.empty()) {
    return raw_result;
  }

  TileEarlyStopResult bound_result = raw_result;
  const std::size_t rows = bound_result.row_passed_flags.size();
  const std::size_t group_size =
      static_cast<std::size_t>(std::max(bind_group_size, 1));

  bound_result.rows_total = rows;
  bound_result.rows_passed = 0;
  bound_result.all_rows_passed = false;

  for (std::size_t group_begin = 0; group_begin < rows; group_begin += group_size) {
    const std::size_t group_end = std::min(group_begin + group_size, rows);
    bool group_passed = true;
    for (std::size_t row = group_begin; row < group_end; ++row) {
      group_passed = group_passed && raw_result.row_passed_flags[row];
    }

    for (std::size_t row = group_begin; row < group_end; ++row) {
      bound_result.row_passed_flags[row] = group_passed;
      if (row < bound_result.row_details.size()) {
        bound_result.row_details[row].row_passed = group_passed;
      }
      if (group_passed) {
        ++bound_result.rows_passed;
      }
    }
  }

  bound_result.all_rows_passed =
      (rows > 0) && (bound_result.rows_passed == rows);
  return bound_result;
}

}  // namespace newcode
