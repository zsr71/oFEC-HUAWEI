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

inline const char* bypass_scheme_name(int scheme_id) {
  switch (scheme_id) {
    case 1:
      return "scheme1";
    case 2:
      return "scheme2";
    default:
      throw std::invalid_argument("bypass_scheme_name: scheme_id must be 1 or 2");
  }
}

inline const std::vector<newcode::mux::MuxEdge>& bypass_edges_for_scheme(
    int scheme_id) {
  switch (scheme_id) {
    case 1:
      return kBypassScheme1Edges;
    case 2:
      return kBypassScheme2Edges;
    default:
      throw std::invalid_argument(
          "bypass_edges_for_scheme: scheme_id must be 1 or 2");
  }
}

}  // namespace app_mux
