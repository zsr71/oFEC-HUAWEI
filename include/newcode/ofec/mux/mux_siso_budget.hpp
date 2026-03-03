#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace newcode::mux {

int pick_siso_active_for_tile(const std::vector<int>& list, std::size_t t);

void apply_siso_budget_g1(std::vector<uint8_t>& state,
                          int siso_active_for_tile);

}  // namespace newcode::mux

