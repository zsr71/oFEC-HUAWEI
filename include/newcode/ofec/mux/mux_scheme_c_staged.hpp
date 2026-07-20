#pragma once

#include <vector>

#include "newcode/ofec/mux/mux_topology.hpp"

namespace newcode::mux {

struct SchemeCResult {
  std::vector<int> stage1_code_to_siso;
  std::vector<int> final_code_to_siso;
  std::vector<int> waiting_codes;
};

SchemeCResult schedule_scheme_c_staged_cpp(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<int>& active_codes,
    const std::vector<int>& free_siso,
    const std::vector<MuxEdge>& extra_bypass_edges);

}  // namespace newcode::mux
