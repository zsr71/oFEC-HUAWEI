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
ValidationResult validate_hiho_active_list(const std::vector<int>& list,
                                           std::size_t tiles_per_win);
ValidationResult validate_siso_active_prefix(const std::vector<int>& list,
                                             std::size_t tiles_to_use);
ValidationResult validate_hiho_active_prefix(const std::vector<int>& list,
                                             std::size_t tiles_to_use);

}  // namespace newcode::mux
