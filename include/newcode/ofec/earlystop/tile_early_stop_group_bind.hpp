#pragma once

#include "newcode/ofec/earlystop/tile_early_stop_result.hpp"

namespace newcode {

TileEarlyStopResult apply_group_bound_early_stop(
    const TileEarlyStopResult& raw_result,
    int bind_group_size);

}  // namespace newcode

#include "ofec/earlystop/tile_early_stop_group_bind.ipp"
