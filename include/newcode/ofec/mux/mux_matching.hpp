#pragma once

#include <cstdint>
#include <vector>

namespace newcode::mux {

struct MatchingResult {
  std::vector<int> code_to_siso;
  std::vector<int> waiting_codes;
};

bool try_augment(int code_idx,
                 const std::vector<std::vector<int>>& adjacency,
                 const std::vector<uint8_t>& free_siso_mask,
                 std::vector<int>& match_siso_to_code,
                 std::vector<uint8_t>& visited_siso);

MatchingResult maximum_bipartite_matching(
    const std::vector<int>& active_codes,
    const std::vector<std::vector<int>>& adjacency,
    const std::vector<int>& free_siso);

MatchingResult augment_from_seed_matching(
    const std::vector<int>& active_codes,
    const std::vector<std::vector<int>>& adjacency,
    const std::vector<int>& free_siso,
    const std::vector<int>& seed_code_to_siso);

}  // namespace newcode::mux

