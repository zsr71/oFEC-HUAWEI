#include "newcode/ofec/mux/mux_scheme_c_staged.hpp"

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace newcode::mux {
namespace {

void validate_scheme_c_args(int n_code, int n_siso, int group_g) {
  if (n_code < 0 || n_siso < 0) {
    throw std::invalid_argument(
        "schedule_scheme_c_staged_cpp: sizes must be >= 0");
  }
  if (group_g < 1) {
    throw std::invalid_argument(
        "schedule_scheme_c_staged_cpp: group_g must be >= 1");
  }
  if (n_code == 0 || n_siso == 0) {
    return;
  }
  if (n_code % group_g != 0 || n_siso % group_g != 0) {
    throw std::invalid_argument(
        "schedule_scheme_c_staged_cpp: n_code and n_siso must be divisible by group_g");
  }
  if (n_code / group_g < 2) {
    throw std::invalid_argument(
        "schedule_scheme_c_staged_cpp: code_per_group must be >= 2");
  }
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

std::vector<uint8_t> build_active_code_mask(int n_code,
                                            const std::vector<int>& active_codes) {
  std::vector<uint8_t> mask(static_cast<std::size_t>(std::max(0, n_code)), 0u);
  for (int code_idx : active_codes) {
    if (code_idx < 0 || code_idx >= n_code) {
      continue;
    }
    mask[static_cast<std::size_t>(code_idx)] = 1u;
  }
  return mask;
}

std::vector<std::vector<int>> build_bypass_map(
    int n_code,
    int n_siso,
    const std::vector<MuxEdge>& extra_bypass_edges) {
  std::vector<std::vector<int>> bypass_map(
      static_cast<std::size_t>(std::max(0, n_code)));
  for (const MuxEdge& edge : extra_bypass_edges) {
    if (edge.code_idx < 0 || edge.code_idx >= n_code) {
      continue;
    }
    if (edge.siso_idx < 0 || edge.siso_idx >= n_siso) {
      continue;
    }
    bypass_map[static_cast<std::size_t>(edge.code_idx)].push_back(edge.siso_idx);
  }
  for (auto& neighbors : bypass_map) {
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
  }
  return bypass_map;
}

std::vector<int> build_local_siso_list_for_group(int group_idx, int siso_per_group) {
  std::vector<int> local_sisos;
  local_sisos.reserve(static_cast<std::size_t>(std::max(0, siso_per_group)));
  const int begin = group_idx * siso_per_group;
  const int end = begin + siso_per_group;
  for (int siso_idx = begin; siso_idx < end; ++siso_idx) {
    local_sisos.push_back(siso_idx);
  }
  return local_sisos;
}

std::vector<int> build_normal_code_list_for_group(int group_idx, int code_per_group) {
  std::vector<int> normal_codes;
  const int normal_count = code_per_group - 2;
  normal_codes.reserve(static_cast<std::size_t>(std::max(0, normal_count)));
  const int begin = group_idx * code_per_group;
  const int end = begin + normal_count;
  for (int code_idx = begin; code_idx < end; ++code_idx) {
    normal_codes.push_back(code_idx);
  }
  return normal_codes;
}

std::vector<int> build_tail_code_list_for_group(int group_idx, int code_per_group) {
  std::vector<int> tail_codes;
  const int tail_count = std::max(2, code_per_group / 4);
  tail_codes.reserve(static_cast<std::size_t>(tail_count));
  const int begin = group_idx * code_per_group + (code_per_group - tail_count);
  const int end = group_idx * code_per_group + code_per_group;
  for (int code_idx = begin; code_idx < end; ++code_idx) {
    tail_codes.push_back(code_idx);
  }
  return tail_codes;
}

bool all_sisos_used(const std::vector<int>& local_sisos,
                    const std::vector<uint8_t>& free_siso_mask,
                    const std::vector<uint8_t>& used_siso) {
  for (int siso_idx : local_sisos) {
    if (siso_idx < 0 || static_cast<std::size_t>(siso_idx) >= used_siso.size()) {
      continue;
    }
    if (free_siso_mask[static_cast<std::size_t>(siso_idx)] &&
        !used_siso[static_cast<std::size_t>(siso_idx)]) {
      return false;
    }
  }
  return true;
}

bool try_assign_first_free(int code_idx,
                           const std::vector<int>& candidate_sisos,
                           const std::vector<uint8_t>& free_siso_mask,
                           std::vector<uint8_t>& used_siso,
                           std::vector<int>& code_to_siso) {
  for (int siso_idx : candidate_sisos) {
    if (siso_idx < 0 ||
        static_cast<std::size_t>(siso_idx) >= free_siso_mask.size()) {
      continue;
    }
    if (!free_siso_mask[static_cast<std::size_t>(siso_idx)] ||
        used_siso[static_cast<std::size_t>(siso_idx)]) {
      continue;
    }
    code_to_siso[static_cast<std::size_t>(code_idx)] = siso_idx;
    used_siso[static_cast<std::size_t>(siso_idx)] = 1u;
    return true;
  }
  return false;
}

std::vector<int> build_waiting_codes(const std::vector<int>& active_codes,
                                     const std::vector<int>& final_code_to_siso) {
  std::vector<int> waiting_codes;
  waiting_codes.reserve(active_codes.size());
  for (int code_idx : active_codes) {
    if (code_idx < 0 ||
        static_cast<std::size_t>(code_idx) >= final_code_to_siso.size()) {
      continue;
    }
    if (final_code_to_siso[static_cast<std::size_t>(code_idx)] < 0) {
      waiting_codes.push_back(code_idx);
    }
  }
  return waiting_codes;
}

}  // namespace

SchemeCResult schedule_scheme_c_staged_cpp(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<int>& active_codes,
    const std::vector<int>& free_siso,
    const std::vector<MuxEdge>& extra_bypass_edges) {
  validate_scheme_c_args(n_code, n_siso, group_g);

  SchemeCResult result;
  result.stage1_code_to_siso.assign(static_cast<std::size_t>(std::max(0, n_code)), -1);
  result.final_code_to_siso.assign(static_cast<std::size_t>(std::max(0, n_code)), -1);
  if (n_code == 0 || n_siso == 0) {
    result.waiting_codes = build_waiting_codes(active_codes, result.final_code_to_siso);
    return result;
  }

  const int code_per_group = n_code / group_g;
  const int siso_per_group = n_siso / group_g;
  const auto free_siso_mask = build_free_siso_mask(n_siso, free_siso);
  const auto active_code_mask = build_active_code_mask(n_code, active_codes);
  const auto bypass_map = build_bypass_map(n_code, n_siso, extra_bypass_edges);
  std::vector<uint8_t> used_siso(static_cast<std::size_t>(n_siso), 0u);

  // Phase 1: normal codes only, local SISO only.
  for (int group_idx = 0; group_idx < group_g; ++group_idx) {
    const auto local_sisos =
        build_local_siso_list_for_group(group_idx, siso_per_group);
    const auto normal_codes =
        build_normal_code_list_for_group(group_idx, code_per_group);
    for (int code_idx : normal_codes) {
      if (!active_code_mask[static_cast<std::size_t>(code_idx)]) {
        continue;
      }
      const bool assigned = try_assign_first_free(
          code_idx, local_sisos, free_siso_mask, used_siso, result.final_code_to_siso);
      if (assigned) {
        result.stage1_code_to_siso[static_cast<std::size_t>(code_idx)] =
            result.final_code_to_siso[static_cast<std::size_t>(code_idx)];
      }
      if (all_sisos_used(local_sisos, free_siso_mask, used_siso)) {
        break;
      }
    }
  }

  // Phase 2: tail codes, local SISO first then bypass edges.
  for (int group_idx = 0; group_idx < group_g; ++group_idx) {
    const auto local_sisos =
        build_local_siso_list_for_group(group_idx, siso_per_group);
    const auto tail_codes =
        build_tail_code_list_for_group(group_idx, code_per_group);
    for (int code_idx : tail_codes) {
      if (!active_code_mask[static_cast<std::size_t>(code_idx)]) {
        continue;
      }
      if (result.final_code_to_siso[static_cast<std::size_t>(code_idx)] >= 0) {
        continue;
      }
      if (try_assign_first_free(
              code_idx, local_sisos, free_siso_mask, used_siso,
              result.final_code_to_siso)) {
        continue;
      }
      (void)try_assign_first_free(
          code_idx, bypass_map[static_cast<std::size_t>(code_idx)], free_siso_mask,
          used_siso, result.final_code_to_siso);
    }
  }

  result.waiting_codes = build_waiting_codes(active_codes, result.final_code_to_siso);
  return result;
}

}  // namespace newcode::mux
