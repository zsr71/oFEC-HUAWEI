#include "newcode/ofec/mux/mux_config_validate.hpp"

#include <sstream>

namespace newcode::mux {

ValidationResult validate_siso_active_list(const std::vector<int>& list,
                                           std::size_t tiles_per_win) {
  if (list.size() != tiles_per_win) {
    std::ostringstream oss;
    oss << "SISO_ACTIVE_LIST length mismatch: got " << list.size()
        << ", expected " << tiles_per_win;
    return ValidationResult{false, oss.str()};
  }

  for (std::size_t i = 0; i < list.size(); ++i) {
    if (list[i] < 0) {
      std::ostringstream oss;
      oss << "SISO_ACTIVE_LIST[" << i << "] must be >= 0, got " << list[i];
      return ValidationResult{false, oss.str()};
    }
  }

  return ValidationResult{};
}

}  // namespace newcode::mux

