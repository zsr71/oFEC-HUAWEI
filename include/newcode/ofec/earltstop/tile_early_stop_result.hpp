#pragma once

#include <cstddef>
#include <vector>

namespace newcode {

struct TileEarlyStopResult {
  bool all_rows_passed = false;
  std::size_t rows_passed = 0;
  std::size_t rows_total = 0;
  std::vector<bool> row_passed_flags;
};

} // namespace newcode
