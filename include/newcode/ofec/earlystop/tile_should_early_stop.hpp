#pragma once

#include "tile_early_stop_stats.hpp"

namespace newcode {

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix);

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix, int detect_mode);

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix,
                            const Params& p,
                            int detect_mode);

} // namespace newcode

#include "ofec/earlystop/tile_should_early_stop.ipp"
