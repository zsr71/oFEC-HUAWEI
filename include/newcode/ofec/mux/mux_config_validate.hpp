#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace newcode::mux {

struct ValidationResult {
  bool ok = true;
  std::string error;
};

ValidationResult validate_siso_active_list(const std::vector<int>& list,
                                           std::size_t tiles_per_win);

}  // namespace newcode::mux

