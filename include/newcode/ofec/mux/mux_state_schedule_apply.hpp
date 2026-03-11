#pragma once

#include <cstdint>
#include <vector>

namespace newcode::mux {

std::vector<int> collect_active_codes_from_state(const std::vector<uint8_t>& state);

std::vector<int> build_free_siso_list(int siso_active_for_tile);

void apply_schedule_result_to_mux_state(std::vector<uint8_t>& state,
                                        const std::vector<int>& final_code_to_siso);

}  // namespace newcode::mux

