#include "newcode/ofec/mux/mux_group_config_validate.hpp"

#include <sstream>

namespace newcode::mux {

ValidationResult validate_mux_group_g(int group_g, std::size_t code_count) {
  if (group_g < 1) {
    std::ostringstream oss;
    oss << "MUX_GROUP_G must be >= 1, got " << group_g;
    return ValidationResult{false, oss.str()};
  }

  if (code_count == 0) {
    return ValidationResult{};
  }

  if (static_cast<std::size_t>(group_g) > code_count) {
    std::ostringstream oss;
    oss << "MUX_GROUP_G must be <= code_count, got " << group_g
        << ", code_count=" << code_count;
    return ValidationResult{false, oss.str()};
  }

  return ValidationResult{};
}

ValidationResult validate_mux_group_runtime(int group_g,
                                            int siso_active_for_tile,
                                            std::size_t code_count) {
  const auto group_ok = validate_mux_group_g(group_g, code_count);
  if (!group_ok.ok) {
    return group_ok;
  }

  if (siso_active_for_tile < 0) {
    std::ostringstream oss;
    oss << "siso_active_for_tile must be >= 0, got " << siso_active_for_tile;
    return ValidationResult{false, oss.str()};
  }

  return ValidationResult{};
}

ValidationResult validate_mux_reconfig_runtime(int group_g,
                                               int siso_active_for_tile,
                                               std::size_t code_count) {
  const auto runtime_ok =
      validate_mux_group_runtime(group_g, siso_active_for_tile, code_count);
  if (!runtime_ok.ok) {
    return runtime_ok;
  }

  if (code_count == 0) {
    return ValidationResult{};
  }

  if (code_count % static_cast<std::size_t>(group_g) != 0) {
    std::ostringstream oss;
    oss << "reconfig mode requires code_count divisible by MUX_GROUP_G, got "
        << "code_count=" << code_count << ", MUX_GROUP_G=" << group_g;
    return ValidationResult{false, oss.str()};
  }

  if (siso_active_for_tile % group_g != 0) {
    std::ostringstream oss;
    oss << "reconfig mode requires siso_active_for_tile divisible by "
           "MUX_GROUP_G, got siso_active_for_tile="
        << siso_active_for_tile << ", MUX_GROUP_G=" << group_g;
    return ValidationResult{false, oss.str()};
  }

  return ValidationResult{};
}

}  // namespace newcode::mux
