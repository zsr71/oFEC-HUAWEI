#include "newcode/ofec/mux/mux_config_validate.hpp"

#include <sstream>

namespace newcode::mux {
namespace {

ValidationResult validate_active_prefix(const std::vector<int>& list,
                                        std::size_t tiles_to_use,
                                        const char* name) {
  if (list.size() < tiles_to_use) {
    std::ostringstream oss;
    oss << name << " length mismatch: got " << list.size()
        << ", expected at least " << tiles_to_use;
    return ValidationResult{false, oss.str()};
  }

  for (std::size_t i = 0; i < tiles_to_use; ++i) {
    if (list[i] < 0) {
      std::ostringstream oss;
      oss << name << "[" << i << "] must be >= 0, got " << list[i];
      return ValidationResult{false, oss.str()};
    }
  }

  return ValidationResult{};
}

ValidationResult validate_active_list(const std::vector<int>& list,
                                      std::size_t tiles_per_win,
                                      const char* name) {
  if (list.size() != tiles_per_win) {
    std::ostringstream oss;
    oss << name << " length mismatch: got " << list.size()
        << ", expected " << tiles_per_win;
    return ValidationResult{false, oss.str()};
  }
  return validate_active_prefix(list, tiles_per_win, name);
}

}  // namespace

ValidationResult validate_siso_active_list(const std::vector<int>& list,
                                           std::size_t tiles_per_win) {
  return validate_active_list(list, tiles_per_win, "SISO_ACTIVE_LIST");
}

ValidationResult validate_hiho_active_list(const std::vector<int>& list,
                                           std::size_t tiles_per_win) {
  return validate_active_list(list, tiles_per_win, "HIHO_ACTIVE_LIST");
}

ValidationResult validate_siso_active_prefix(const std::vector<int>& list,
                                             std::size_t tiles_to_use) {
  return validate_active_prefix(list, tiles_to_use, "SISO_ACTIVE_LIST");
}

ValidationResult validate_hiho_active_prefix(const std::vector<int>& list,
                                             std::size_t tiles_to_use) {
  return validate_active_prefix(list, tiles_to_use, "HIHO_ACTIVE_LIST");
}

}  // namespace newcode::mux
