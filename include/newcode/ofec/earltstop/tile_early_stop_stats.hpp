#pragma once

#include "newcode/llr_utils.hpp"
#include "newcode/common/matrix/matrix.hpp"

#include "tile_early_stop_result.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace newcode {

template <typename LLR>
TileEarlyStopResult tile_early_stop_stats(const matrix::Matrix<LLR>& lin_matrix);

} // namespace newcode

#include "ofec/earltstop/tile_early_stop_stats.ipp"
