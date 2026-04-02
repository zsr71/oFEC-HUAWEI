#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace newcode {

struct TileEarlyStopRowDetail {
  bool row_passed = false;
  bool bch_passed = false;
  bool overall_parity_passed = false;
  std::uint8_t syndrome_bits = 0;
};

struct TileEarlyStopResult {
  bool all_rows_passed = false;
  std::size_t rows_passed = 0;
  std::size_t rows_total = 0;
  std::vector<bool> row_passed_flags;
  std::vector<TileEarlyStopRowDetail> row_details;
};

} // namespace newcode
