#include "newcode/ofec/mux/mux_matching.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace newcode::mux {
namespace {

int infer_siso_count(const std::vector<std::vector<int>>& adjacency,
                     const std::vector<int>& free_siso) {
  int max_siso = -1;
  for (const auto& neighbors : adjacency) {
    for (int siso_idx : neighbors) {
      max_siso = std::max(max_siso, siso_idx);
    }
  }
  for (int siso_idx : free_siso) {
    max_siso = std::max(max_siso, siso_idx);
  }
  return max_siso + 1;
}

std::vector<uint8_t> build_free_siso_mask(int n_siso,
                                          const std::vector<int>& free_siso) {
  std::vector<uint8_t> mask(static_cast<std::size_t>(std::max(0, n_siso)), 0u);
  for (int siso_idx : free_siso) {
    if (siso_idx < 0 || siso_idx >= n_siso) {
      continue;
    }
    mask[static_cast<std::size_t>(siso_idx)] = 1u;
  }
  return mask;
}

MatchingResult build_matching_result(const std::vector<int>& active_codes,
                                     const std::vector<int>& match_siso_to_code,
                                     std::size_t n_code) {
  MatchingResult result;
  result.code_to_siso.assign(n_code, -1);
  for (std::size_t siso_idx = 0; siso_idx < match_siso_to_code.size(); ++siso_idx) {
    const int code_idx = match_siso_to_code[siso_idx];
    if (code_idx < 0 || static_cast<std::size_t>(code_idx) >= n_code) {
      continue;
    }
    result.code_to_siso[static_cast<std::size_t>(code_idx)] =
        static_cast<int>(siso_idx);
  }
  for (int code_idx : active_codes) {
    if (code_idx < 0 || static_cast<std::size_t>(code_idx) >= n_code) {
      continue;
    }
    if (result.code_to_siso[static_cast<std::size_t>(code_idx)] < 0) {
      result.waiting_codes.push_back(code_idx);
    }
  }
  return result;
}

}  // namespace

bool try_augment(int code_idx,
                 const std::vector<std::vector<int>>& adjacency,
                 const std::vector<uint8_t>& free_siso_mask,
                 std::vector<int>& match_siso_to_code,
                 std::vector<uint8_t>& visited_siso) {
  if (code_idx < 0 || static_cast<std::size_t>(code_idx) >= adjacency.size()) {
    return false;
  }
  for (int siso_idx : adjacency[static_cast<std::size_t>(code_idx)]) {
    if (siso_idx < 0 || static_cast<std::size_t>(siso_idx) >= free_siso_mask.size()) {
      continue;
    }
    if (!free_siso_mask[static_cast<std::size_t>(siso_idx)] ||
        visited_siso[static_cast<std::size_t>(siso_idx)]) {
      continue;
    }

    visited_siso[static_cast<std::size_t>(siso_idx)] = 1u;
    const int holder = match_siso_to_code[static_cast<std::size_t>(siso_idx)];
    if (holder < 0 ||
        try_augment(holder, adjacency, free_siso_mask, match_siso_to_code, visited_siso)) {
      match_siso_to_code[static_cast<std::size_t>(siso_idx)] = code_idx;
      return true;
    }
  }
  return false;
}

MatchingResult maximum_bipartite_matching(
    const std::vector<int>& active_codes,
    const std::vector<std::vector<int>>& adjacency,
    const std::vector<int>& free_siso) {
  const int n_siso = infer_siso_count(adjacency, free_siso);
  const auto free_siso_mask = build_free_siso_mask(n_siso, free_siso);
  std::vector<int> match_siso_to_code(
      static_cast<std::size_t>(std::max(0, n_siso)), -1);

  for (int code_idx : active_codes) {
    std::vector<uint8_t> visited_siso(
        static_cast<std::size_t>(std::max(0, n_siso)), 0u);
    try_augment(code_idx, adjacency, free_siso_mask, match_siso_to_code, visited_siso);
  }

  return build_matching_result(active_codes, match_siso_to_code, adjacency.size());
}

MatchingResult augment_from_seed_matching(
    const std::vector<int>& active_codes,
    const std::vector<std::vector<int>>& adjacency,
    const std::vector<int>& free_siso,
    const std::vector<int>& seed_code_to_siso) {
  const int n_siso = infer_siso_count(adjacency, free_siso);
  const auto free_siso_mask = build_free_siso_mask(n_siso, free_siso);
  std::vector<int> match_siso_to_code(
      static_cast<std::size_t>(std::max(0, n_siso)), -1);

  for (std::size_t code_idx = 0; code_idx < seed_code_to_siso.size(); ++code_idx) {
    const int siso_idx = seed_code_to_siso[code_idx];
    if (siso_idx < 0 || siso_idx >= n_siso) {
      continue;
    }
    if (!free_siso_mask[static_cast<std::size_t>(siso_idx)]) {
      continue;
    }
    match_siso_to_code[static_cast<std::size_t>(siso_idx)] =
        static_cast<int>(code_idx);
  }

  for (int code_idx : active_codes) {
    if (code_idx < 0 || static_cast<std::size_t>(code_idx) >= seed_code_to_siso.size()) {
      continue;
    }
    if (seed_code_to_siso[static_cast<std::size_t>(code_idx)] >= 0) {
      continue;
    }
    std::vector<uint8_t> visited_siso(
        static_cast<std::size_t>(std::max(0, n_siso)), 0u);
    try_augment(code_idx, adjacency, free_siso_mask, match_siso_to_code, visited_siso);
  }

  return build_matching_result(active_codes, match_siso_to_code, adjacency.size());
}

}  // namespace newcode::mux

