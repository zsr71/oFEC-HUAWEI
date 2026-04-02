#include "newcode/ofec/mux/mux_group_budget.hpp"

#include "newcode/ofec/mux/mux_siso_budget.hpp"
#include "newcode/ofec/mux/mux_state_builder.hpp"

#include <algorithm>
#include <bit>
#include <stdexcept>
#include <vector>

namespace newcode::mux {
namespace {

void validate_grouped_budget_args(std::size_t code_count,
                                  int siso_active_for_tile,
                                  int group_g) {
  if (group_g < 1) {
    throw std::invalid_argument(
        "apply_siso_budget_grouped: group_g must be >= 1");
  }
  if (siso_active_for_tile < 0) {
    throw std::invalid_argument(
        "apply_siso_budget_grouped: siso_active_for_tile must be >= 0");
  }
  if (code_count == 0) {
    return;
  }
  if (static_cast<std::size_t>(group_g) > code_count) {
    throw std::invalid_argument(
        "apply_siso_budget_grouped: group_g must be <= state.size()");
  }
}

std::vector<std::size_t> collect_need_indices_in_range(
    const std::vector<uint8_t>& state,
    std::size_t begin,
    std::size_t end) {
  std::vector<std::size_t> need_indices;
  need_indices.reserve(end - begin);
  for (std::size_t i = begin; i < end; ++i) {
    if (state[i] == static_cast<uint8_t>(StateTag::NeedSiso)) {
      need_indices.push_back(i);
    }
  }
  return need_indices;
}

void trim_need_indices_by_budget(std::vector<uint8_t>& state,
                                 const std::vector<std::size_t>& need_indices,
                                 int keep_budget) {
  if (keep_budget < 0) {
    throw std::invalid_argument(
        "trim_need_indices_by_budget: keep_budget must be >= 0");
  }

  const std::size_t start = std::min<std::size_t>(
      static_cast<std::size_t>(keep_budget), need_indices.size());
  for (std::size_t k = start; k < need_indices.size(); ++k) {
    state[need_indices[k]] = static_cast<uint8_t>(StateTag::Unscheduled);
  }
}

int syndrome_nonzero_count(std::uint8_t syndrome_bits) {
  return std::popcount(static_cast<unsigned int>(syndrome_bits));
}

int early_stop_priority_score(const TileEarlyStopRowDetail& detail,
                              int priority_rule) {
  const int syndrome_count = syndrome_nonzero_count(detail.syndrome_bits);
  switch (priority_rule) {
    case 0:  // harder_first
      return syndrome_count * 10 + (detail.overall_parity_passed ? 0 : 1);
    case 1:  // near_threshold_first
      if (detail.bch_passed && !detail.overall_parity_passed) {
        return 100;
      }
      return (detail.bch_passed ? 50 : 0) +
             (detail.overall_parity_passed ? 5 : 0) - syndrome_count;
    default:
      return 0;
  }
}

void trim_need_indices_by_priority(std::vector<uint8_t>& state,
                                   std::vector<std::size_t> need_indices,
                                   int keep_budget,
                                   const TileEarlyStopResult& stats,
                                   int priority_rule) {
  if (keep_budget < 0) {
    throw std::invalid_argument(
        "trim_need_indices_by_priority: keep_budget must be >= 0");
  }

  std::stable_sort(
      need_indices.begin(),
      need_indices.end(),
      [&](std::size_t lhs, std::size_t rhs) {
        const auto lhs_score =
            early_stop_priority_score(stats.row_details[lhs], priority_rule);
        const auto rhs_score =
            early_stop_priority_score(stats.row_details[rhs], priority_rule);
        return lhs_score > rhs_score;
      });

  const std::size_t start = std::min<std::size_t>(
      static_cast<std::size_t>(keep_budget), need_indices.size());
  for (std::size_t k = start; k < need_indices.size(); ++k) {
    state[need_indices[k]] = static_cast<uint8_t>(StateTag::Unscheduled);
  }
}

}  // namespace

std::vector<GroupRange> build_groups_even(std::size_t total, int group_g) {
  if (group_g < 1) {
    throw std::invalid_argument("build_groups_even: group_g must be >= 1");
  }
  if (total > 0 && static_cast<std::size_t>(group_g) > total) {
    throw std::invalid_argument("build_groups_even: group_g must be <= total");
  }

  const std::size_t g = static_cast<std::size_t>(group_g);
  const std::size_t base = total / g;
  const std::size_t rem = total % g;

  std::vector<GroupRange> groups;
  groups.reserve(g);
  std::size_t start = 0;
  for (std::size_t i = 0; i < g; ++i) {
    const std::size_t span = base + (i < rem ? 1u : 0u);
    groups.push_back(GroupRange{start, start + span});
    start += span;
  }
  return groups;
}

std::vector<int> split_budget_even(int total_budget, int group_g) {
  if (total_budget < 0) {
    throw std::invalid_argument("split_budget_even: total_budget must be >= 0");
  }
  if (group_g < 1) {
    throw std::invalid_argument("split_budget_even: group_g must be >= 1");
  }

  const int base = total_budget / group_g;
  const int rem = total_budget % group_g;

  std::vector<int> budgets(static_cast<std::size_t>(group_g), base);
  for (int i = 0; i < rem; ++i) {
    budgets[static_cast<std::size_t>(i)] += 1;
  }
  return budgets;
}

void apply_siso_budget_grouped(std::vector<uint8_t>& state,
                               int siso_active_for_tile,
                               int group_g) {
  if (state.empty()) {
    return;
  }

  validate_grouped_budget_args(state.size(), siso_active_for_tile, group_g);
  if (group_g <= 1) {
    apply_siso_budget_g1(state, siso_active_for_tile);
    return;
  }

  const auto groups = build_groups_even(state.size(), group_g);
  const auto budgets = split_budget_even(siso_active_for_tile, group_g);
  for (int g = 0; g < group_g; ++g) {
    const GroupRange range = groups[static_cast<std::size_t>(g)];
    const auto need_indices =
        collect_need_indices_in_range(state, range.begin, range.end);
    trim_need_indices_by_budget(
        state, need_indices, budgets[static_cast<std::size_t>(g)]);
  }
}

void apply_siso_budget_grouped_priority(std::vector<uint8_t>& state,
                                        int siso_active_for_tile,
                                        int group_g,
                                        const TileEarlyStopResult& stats,
                                        int priority_rule) {
  if (state.empty()) {
    return;
  }

  validate_grouped_budget_args(state.size(), siso_active_for_tile, group_g);
  if (stats.row_details.size() != state.size()) {
    apply_siso_budget_grouped(state, siso_active_for_tile, group_g);
    return;
  }
  if (group_g <= 1) {
    apply_siso_budget_g1_priority(state,
                                  siso_active_for_tile,
                                  stats,
                                  priority_rule);
    return;
  }

  const auto groups = build_groups_even(state.size(), group_g);
  const auto budgets = split_budget_even(siso_active_for_tile, group_g);
  for (int g = 0; g < group_g; ++g) {
    const GroupRange range = groups[static_cast<std::size_t>(g)];
    auto need_indices =
        collect_need_indices_in_range(state, range.begin, range.end);
    trim_need_indices_by_priority(state,
                                  std::move(need_indices),
                                  budgets[static_cast<std::size_t>(g)],
                                  stats,
                                  priority_rule);
  }
}

}  // namespace newcode::mux
