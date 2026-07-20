#pragma once

#include <stdexcept>
#include <vector>

#include "newcode/ofec/mux/mux_topology.hpp"

namespace app_mux {

inline const std::vector<newcode::mux::MuxEdge> kBypassScheme1Edges = {
    {6, 11},  {6, 15},  {7, 11},  {7, 15},
    {14, 11}, {14, 15}, {15, 11}, {15, 15},
    {22, 3},  {22, 7},  {23, 3},  {23, 7},
    {30, 3},  {30, 7},  {31, 3},  {31, 7},
};

inline const std::vector<newcode::mux::MuxEdge> kBypassScheme2Edges = {
    {6, 11},  {6, 15},  {7, 11},  {7, 15},
    {14, 11}, {14, 15}, {15, 11}, {15, 15},
    {22, 3},  {22, 7},  {23, 3},  {23, 7},
    {30, 3},  {30, 7},  {31, 3},  {31, 7},
    {6, 7},   {7, 7},   {14, 3},  {15, 3},
    {22, 15}, {23, 15}, {30, 11}, {31, 11},
};

// scheme3 口径：
// - 整体按 2 组理解，每组 16 个 code
// - 每组本地可用 8 个 SISO
// - 每组最后 4 个 code 允许旁路到另一组的后 4/前 4 个 SISO
// 该方案在 group_g=2 且 siso_active_for_tile=16 时最符合设计意图。
inline const std::vector<newcode::mux::MuxEdge> kBypassScheme3Edges = {
    {12, 12}, {12, 13}, {12, 14}, {12, 15},
    {13, 12}, {13, 13}, {13, 14}, {13, 15},
    {14, 12}, {14, 13}, {14, 14}, {14, 15},
    {15, 12}, {15, 13}, {15, 14}, {15, 15},
    {28, 4},  {28, 5},  {28, 6},  {28, 7},
    {29, 4},  {29, 5},  {29, 6},  {29, 7},
    {30, 4},  {30, 5},  {30, 6},  {30, 7},
    {31, 4},  {31, 5},  {31, 6},  {31, 7},
};

inline const char* bypass_scheme_name(int scheme_id) {
  switch (scheme_id) {
    case 1:
      return "scheme1";
    case 2:
      return "scheme2";
    case 3:
      return "scheme3";
    default:
      throw std::invalid_argument("bypass_scheme_name: scheme_id must be 1, 2 or 3");
  }
}

inline const std::vector<newcode::mux::MuxEdge>& bypass_edges_for_scheme(
    int scheme_id) {
  switch (scheme_id) {
    case 1:
      return kBypassScheme1Edges;
    case 2:
      return kBypassScheme2Edges;
    case 3:
      return kBypassScheme3Edges;
    default:
      throw std::invalid_argument(
          "bypass_edges_for_scheme: scheme_id must be 1, 2 or 3");
  }
}

}  // namespace app_mux
