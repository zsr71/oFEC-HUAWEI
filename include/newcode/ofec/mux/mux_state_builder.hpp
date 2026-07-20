#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "newcode/ofec/earlystop/tile_early_stop_result.hpp"

namespace newcode::mux {

enum class StateTag : uint8_t {
  NeedSiso = 0,
  EarlyStopped = 1,
  Unscheduled = 2
};

std::vector<uint8_t> build_state_from_early_stop(
    const TileEarlyStopResult& stats);

struct StateCount {
  std::size_t n0_need_siso = 0;
  std::size_t n1_early_stopped = 0;
  std::size_t n2_unscheduled = 0;
};

StateCount count_state012(const std::vector<uint8_t>& state);

}  // namespace newcode::mux

