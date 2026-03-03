#include "newcode/ofec/mux/mux_siso_budget.hpp"
#include "newcode/ofec/mux/mux_state_builder.hpp"

#include <stdexcept>
#include <vector>

namespace newcode::mux {

int pick_siso_active_for_tile(const std::vector<int>& list, std::size_t t) {
  if (t >= list.size()) {
    throw std::out_of_range("pick_siso_active_for_tile: tile index out of range");
  }
  return list[t];
}

void apply_siso_budget_g1(std::vector<uint8_t>& state,
                          int siso_active_for_tile) {
  if (siso_active_for_tile < 0) {
    throw std::invalid_argument("apply_siso_budget_g1: siso_active_for_tile must be >= 0");
  }

  std::vector<std::size_t> need_indices;
  need_indices.reserve(state.size());

  for (std::size_t i = 0; i < state.size(); ++i) {
    if (state[i] == static_cast<uint8_t>(StateTag::NeedSiso)) {
      need_indices.push_back(i);
    }
  }

  if (static_cast<int>(need_indices.size()) <= siso_active_for_tile) {
    return;
  }

  for (std::size_t k = static_cast<std::size_t>(siso_active_for_tile);
       k < need_indices.size();
       ++k) {
    state[need_indices[k]] = static_cast<uint8_t>(StateTag::Unscheduled);
  }
}

}  // namespace newcode::mux
