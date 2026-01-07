#pragma once

#include "tile_early_stop_stats.hpp"

namespace newcode {

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix);

} // namespace newcode

#include "ofec/earlystop/tile_should_early_stop.ipp"
