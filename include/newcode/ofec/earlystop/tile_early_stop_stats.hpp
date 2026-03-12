#pragma once

#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/common/matrix/matrix.hpp"

#include "tile_early_stop_result.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>

namespace newcode {

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop_v1(const matrix::Matrix<LLR>& lin_matrix);

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop_v2(const matrix::Matrix<LLR>& lin_matrix);

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop(const matrix::Matrix<LLR>& lin_matrix,
                                           int detect_mode) {
  switch (detect_mode) {
    case 1:
      return detect_tile_early_stop_v1(lin_matrix);
    case 2:
      return detect_tile_early_stop_v2(lin_matrix);
    default:
      throw std::invalid_argument(
          "detect_tile_early_stop: detect_mode must be 1 or 2");
  }
}

} // namespace newcode

#include "ofec/earlystop/tile_early_stop_stats.ipp"
#include "ofec/earlystop/tile_early_stop_stats2.ipp"
