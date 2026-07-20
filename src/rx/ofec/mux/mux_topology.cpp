#include "newcode/ofec/mux/mux_topology.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>
#include <vector>

namespace newcode::mux {
namespace {

struct GroupShape {
  int code_per_group = 0;
  int siso_per_group = 0;
};

GroupShape validate_group_shape(int n_code, int n_siso, int group_g) {
  if (n_code < 0 || n_siso < 0) {
    throw std::invalid_argument("validate_group_shape: sizes must be >= 0");
  }
  if (group_g < 1) {
    throw std::invalid_argument("validate_group_shape: group_g must be >= 1");
  }
  if (n_code == 0 || n_siso == 0) {
    return GroupShape{};
  }
  if (n_code % group_g != 0 || n_siso % group_g != 0) {
    throw std::invalid_argument(
        "validate_group_shape: n_code and n_siso must be divisible by group_g");
  }
  return GroupShape{n_code / group_g, n_siso / group_g};
}

}  // namespace

std::vector<MuxEdge> build_local_edges(int n_code, int n_siso, int group_g) {
  const GroupShape shape = validate_group_shape(n_code, n_siso, group_g);
  std::vector<MuxEdge> edges;
  if (n_code == 0 || n_siso == 0) {
    return edges;
  }

  edges.reserve(static_cast<std::size_t>(n_code * shape.siso_per_group));
  for (int code_idx = 0; code_idx < n_code; ++code_idx) {
    const int g = code_idx / shape.code_per_group;
    const int siso_begin = g * shape.siso_per_group;
    const int siso_end = siso_begin + shape.siso_per_group;
    for (int siso_idx = siso_begin; siso_idx < siso_end; ++siso_idx) {
      edges.push_back(MuxEdge{code_idx, siso_idx});
    }
  }
  return edges;
}

std::vector<MuxEdge> build_allowed_edges(
    int n_code,
    int n_siso,
    int group_g,
    const std::vector<MuxEdge>& extra_bypass_edges) {
  std::vector<MuxEdge> local_edges = build_local_edges(n_code, n_siso, group_g);
  std::set<std::pair<int, int>> unique_edges;
  for (const MuxEdge& edge : local_edges) {
    unique_edges.emplace(edge.code_idx, edge.siso_idx);
  }
  for (const MuxEdge& edge : extra_bypass_edges) {
    if (edge.code_idx < 0 || edge.siso_idx < 0) {
      continue;
    }
    if (edge.code_idx >= n_code || edge.siso_idx >= n_siso) {
      continue;
    }
    unique_edges.emplace(edge.code_idx, edge.siso_idx);
  }

  std::vector<MuxEdge> edges;
  edges.reserve(unique_edges.size());
  for (const auto& entry : unique_edges) {
    edges.push_back(MuxEdge{entry.first, entry.second});
  }
  return edges;
}

std::vector<std::vector<int>> build_adjacency(const std::vector<MuxEdge>& edges,
                                              int n_code) {
  if (n_code < 0) {
    throw std::invalid_argument("build_adjacency: n_code must be >= 0");
  }
  std::vector<std::vector<int>> adjacency(static_cast<std::size_t>(n_code));
  for (const MuxEdge& edge : edges) {
    if (edge.code_idx < 0 || edge.code_idx >= n_code || edge.siso_idx < 0) {
      continue;
    }
    adjacency[static_cast<std::size_t>(edge.code_idx)].push_back(edge.siso_idx);
  }
  for (auto& neighbors : adjacency) {
    std::sort(neighbors.begin(), neighbors.end());
    neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
  }
  return adjacency;
}

}  // namespace newcode::mux

