#pragma once

#include <stdexcept>

#include "newcode/ofec/earlystop/row_early_stop_process_1.hpp"
#include "newcode/ofec/earlystop/row_early_stop_process_2.hpp"
#include "newcode/ofec/earlystop/row_early_stop_process_3.hpp"

namespace newcode {

template <typename LLR>
bool apply_row_early_stop_action(const LLR* lin256,
                                 const LLR* lch256,
                                 float* y2_256,
                                 const newcode::Params& p) {
  switch (p.EARLY_STOP_ACTION_MODE) {
    case 1:
      row_early_stop_process_1(lin256, lch256, y2_256, p);
      return true;
    case 2:
      row_early_stop_process_2(lin256, lch256, y2_256, p);
      return true;
    case 3:
      return row_early_stop_process_3(lin256, lch256, y2_256, p);
    default:
      throw std::invalid_argument(
          "EARLY_STOP_ACTION_MODE must be 1, 2 or 3 in current implementation");
  }
}

}  // namespace newcode
