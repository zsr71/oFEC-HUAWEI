#pragma once

#include <cstddef>
#include <vector>

#include "newcode/ofec/mux/mux_config_validate.hpp"

namespace newcode::mux {

ValidationResult validate_mux_group_g(int group_g, std::size_t code_count);

int pick_mux_group_g_for_tile(const std::vector<int>& group_g_list,
                              std::size_t tile_index,
                              int fallback_group_g);

ValidationResult validate_mux_group_runtime(int group_g,
                                            int siso_active_for_tile,
                                            std::size_t code_count);

ValidationResult validate_mux_reconfig_runtime(int group_g,
                                               int siso_active_for_tile,
                                               std::size_t code_count);

}  // namespace newcode::mux
