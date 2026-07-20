#include "newcode/ofec/mux/mux_state_schedule_apply.hpp"

#include <cstdint>
#include <stdexcept>
#include <vector>

#include "newcode/ofec/mux/mux_state_builder.hpp"

namespace newcode::mux {

std::vector<int> collect_active_codes_from_state(const std::vector<uint8_t>& state) {
  std::vector<int> active_codes;
  active_codes.reserve(state.size());
  for (std::size_t idx = 0; idx < state.size(); ++idx) {
    if (state[idx] == static_cast<uint8_t>(StateTag::NeedSiso)) {
      active_codes.push_back(static_cast<int>(idx));
    }
  }
  return active_codes;
}

std::vector<int> build_free_siso_list(int siso_active_for_tile) {
  if (siso_active_for_tile < 0) {
    throw std::invalid_argument(
        "build_free_siso_list: siso_active_for_tile must be >= 0");
  }
  std::vector<int> free_siso;
  free_siso.reserve(static_cast<std::size_t>(siso_active_for_tile));
  for (int siso_idx = 0; siso_idx < siso_active_for_tile; ++siso_idx) {
    free_siso.push_back(siso_idx);
  }
  return free_siso;
}

void apply_schedule_result_to_mux_state(std::vector<uint8_t>& state,
                                        const std::vector<int>& final_code_to_siso) {
  for (std::size_t idx = 0; idx < state.size(); ++idx) {
    if (state[idx] != static_cast<uint8_t>(StateTag::NeedSiso)) {
      continue;
    }
    const bool matched =
        idx < final_code_to_siso.size() && final_code_to_siso[idx] >= 0;
    if (!matched) {
      state[idx] = static_cast<uint8_t>(StateTag::Unscheduled);
    }
  }
}

}  // namespace newcode::mux

