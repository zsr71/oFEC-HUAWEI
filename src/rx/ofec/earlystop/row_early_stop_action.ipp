#pragma once

#include <stdexcept>

#include "newcode/ofec/earlystop/row_early_stop_process_1.hpp"
#include "newcode/ofec/earlystop/row_early_stop_process_2.hpp"

namespace newcode {

template <typename LLR>
void apply_row_early_stop_action(const LLR* lin256,
                                 const LLR* lch256,
                                 float* y2_256,
                                 const newcode::Params& p) {
  switch (p.EARLY_STOP_ACTION_MODE) {
    case 1:
      row_early_stop_process_1(lin256, lch256, y2_256, p);
      return;
    case 2:
      row_early_stop_process_2(lin256, lch256, y2_256, p);
      return;
    default:
      throw std::invalid_argument(
          "EARLY_STOP_ACTION_MODE must be 1 or 2 in current implementation");
  }
}

}  // namespace newcode
