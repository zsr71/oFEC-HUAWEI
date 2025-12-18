#pragma once

#include "newcode/bch_255_239.hpp"
#include "newcode/llr_utils.hpp"
#include "newcode/matrix.hpp"

#include "tile_early_stop_result.hpp"

#include <array>
#include <cstddef>
#include <cstdint>

namespace newcode {

template <typename LLR>
TileEarlyStopResult tile_early_stop_stats(const Matrix<LLR>& lin_matrix);

} // namespace newcode

#include "ofec/earltstop/tile_early_stop_stats.ipp"
