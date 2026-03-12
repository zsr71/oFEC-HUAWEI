#pragma once

#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/common/matrix/matrix.hpp"
#include "newcode/params.hpp"

#include "tile_early_stop_result.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace newcode {

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop_v1(const matrix::Matrix<LLR>& lin_matrix,
                                              const Params& p);

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop_v2(const matrix::Matrix<LLR>& lin_matrix,
                                              const Params& p);

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop(const matrix::Matrix<LLR>& lin_matrix,
                                           const Params& p,
                                           int detect_mode) {
  switch (detect_mode) {
    case 1:
      return detect_tile_early_stop_v1(lin_matrix, p);
    case 2:
      return detect_tile_early_stop_v2(lin_matrix, p);
    default:
      throw std::invalid_argument(
          "detect_tile_early_stop: detect_mode must be 1 or 2");
  }
}

} // namespace newcode

#include "ofec/earlystop/tile_early_stop_stats.ipp"
#include "ofec/earlystop/tile_early_stop_stats2.ipp"
