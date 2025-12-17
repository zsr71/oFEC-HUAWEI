#pragma once

#include "ofec_tile_input.ipp"

namespace newcode {
namespace detail {

template <typename LLR>
TileEarlyStopResult run_tile_early_stop(const TilePrepared<LLR>& prep) {
  return tile_early_stop_stats(prep.lin_matrix);
}

} // namespace detail
} // namespace newcode
