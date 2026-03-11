#pragma once

#include <cstddef>

#include "newcode/ofec/mux/mux_config_validate.hpp"

namespace newcode::mux {

ValidationResult validate_mux_group_g(int group_g, std::size_t code_count);

ValidationResult validate_mux_group_runtime(int group_g,
                                            int siso_active_for_tile,
                                            std::size_t code_count);

ValidationResult validate_mux_reconfig_runtime(int group_g,
                                               int siso_active_for_tile,
                                               std::size_t code_count);

}  // namespace newcode::mux
