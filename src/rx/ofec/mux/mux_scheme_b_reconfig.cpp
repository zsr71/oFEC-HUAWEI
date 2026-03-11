#include "newcode/ofec/mux/mux_scheme_b_reconfig.hpp"

#include "newcode/ofec/mux/mux_matching.hpp"
#include "newcode/ofec/mux/mux_topology.hpp"

#include <stdexcept>
#include <vector>

namespace newcode::mux {

SchemeBResult schedule_scheme_b_reconfig_cpp(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<int>& active_codes,
    const std::vector<int>& free_siso,
    const std::vector<MuxEdge>& extra_bypass_edges) {
  if (n_code < 0 || n_siso < 0) {
    throw std::invalid_argument(
        "schedule_scheme_b_reconfig_cpp: sizes must be >= 0");
  }

  const auto local_edges = build_local_edges(n_code, n_siso, group_g);
  const auto all_edges =
      build_allowed_edges(n_code, n_siso, group_g, extra_bypass_edges);
  const auto local_adjacency = build_adjacency(local_edges, n_code);
  const auto all_adjacency = build_adjacency(all_edges, n_code);

  const MatchingResult stage1 =
      maximum_bipartite_matching(active_codes, local_adjacency, free_siso);
  const MatchingResult final = augment_from_seed_matching(
      active_codes, all_adjacency, free_siso, stage1.code_to_siso);

  SchemeBResult result;
  result.stage1_code_to_siso = stage1.code_to_siso;
  result.final_code_to_siso = final.code_to_siso;
  result.waiting_codes = final.waiting_codes;
  return result;
}

}  // namespace newcode::mux

