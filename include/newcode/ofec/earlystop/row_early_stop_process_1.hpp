#pragma once

#include "newcode/params.hpp"

namespace newcode {

template <typename LLR>
void row_early_stop_process_1(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p);

} // namespace newcode

#include "ofec/earlystop/row_early_stop_process_1.ipp"
