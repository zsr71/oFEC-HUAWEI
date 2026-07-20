#include "newcode/ofec/mux/mux_siso_budget.hpp"
#include "newcode/ofec/mux/mux_state_builder.hpp"

#include <algorithm>
#include <bit>
#include <stdexcept>
#include <vector>

namespace newcode::mux {
namespace {

int syndrome_nonzero_count(std::uint8_t syndrome_bits) {
  return std::popcount(static_cast<unsigned int>(syndrome_bits));
}

int early_stop_priority_score(const TileEarlyStopRowDetail& detail,
                              int priority_rule) {
  const int syndrome_count = syndrome_nonzero_count(detail.syndrome_bits);
  switch (priority_rule) {
    case 0: {  // harder_first
      return syndrome_count * 10 + (detail.overall_parity_passed ? 0 : 1);
    }
    case 1: {  // near_threshold_first
      if (detail.bch_passed && !detail.overall_parity_passed) {
        return 100;
      }
      return (detail.bch_passed ? 50 : 0) +
             (detail.overall_parity_passed ? 5 : 0) - syndrome_count;
    }
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
  if (need_indices.empty()) {
    return;
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

  const std::size_t start =
      std::min<std::size_t>(static_cast<std::size_t>(keep_budget),
                            need_indices.size());
  for (std::size_t k = start; k < need_indices.size(); ++k) {
    state[need_indices[k]] = static_cast<uint8_t>(StateTag::Unscheduled);
  }
}

}  // namespace

int pick_siso_active_for_tile(const std::vector<int>& list, std::size_t t) {
  if (t >= list.size()) {
    throw std::out_of_range("pick_siso_active_for_tile: tile index out of range");
  }
  return list[t];
}

int pick_hiho_active_for_tile(const std::vector<int>& list, std::size_t t) {
  if (t >= list.size()) {
    throw std::out_of_range("pick_hiho_active_for_tile: tile index out of range");
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

void apply_siso_budget_g1_priority(std::vector<uint8_t>& state,
                                   int siso_active_for_tile,
                                   const TileEarlyStopResult& stats,
                                   int priority_rule) {
  if (siso_active_for_tile < 0) {
    throw std::invalid_argument(
        "apply_siso_budget_g1_priority: siso_active_for_tile must be >= 0");
  }
  if (stats.row_details.size() != state.size()) {
    apply_siso_budget_g1(state, siso_active_for_tile);
    return;
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

  trim_need_indices_by_priority(
      state, std::move(need_indices), siso_active_for_tile, stats, priority_rule);
}

}  // namespace newcode::mux
