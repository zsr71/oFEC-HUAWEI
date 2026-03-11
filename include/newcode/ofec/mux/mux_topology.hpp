#pragma once

#include <vector>

namespace newcode::mux {

struct MuxEdge {
  int code_idx = -1;
  int siso_idx = -1;
};

std::vector<MuxEdge> build_local_edges(int n_code, int n_siso, int group_g);

std::vector<MuxEdge> build_allowed_edges(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<MuxEdge>& extra_bypass_edges);

std::vector<std::vector<int>> build_adjacency(const std::vector<MuxEdge>& edges,
                                              int n_code);

}  // namespace newcode::mux

