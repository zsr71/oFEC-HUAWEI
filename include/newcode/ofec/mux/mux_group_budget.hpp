#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "newcode/ofec/earlystop/tile_early_stop_result.hpp"

namespace newcode::mux {

struct GroupRange {
  std::size_t begin = 0;  // inclusive
  std::size_t end = 0;    // exclusive
};

std::vector<GroupRange> build_groups_even(std::size_t total, int group_g);

std::vector<int> split_budget_even(int total_budget, int group_g);

void apply_siso_budget_grouped(std::vector<uint8_t>& state,
                               int siso_active_for_tile,
                               int group_g);

void apply_siso_budget_grouped_priority(std::vector<uint8_t>& state,
                                        int siso_active_for_tile,
                                        int group_g,
                                        const TileEarlyStopResult& stats,
                                        int priority_rule);

}  // namespace newcode::mux
