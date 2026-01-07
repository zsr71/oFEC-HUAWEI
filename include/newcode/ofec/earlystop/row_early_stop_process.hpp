#pragma once

#include "newcode/params.hpp"

namespace newcode {

template <typename LLR>
void row_early_stop_process(const LLR* lin256,
                            float* y2_256,
                            const newcode::Params& p);

} // namespace newcode

#include "ofec/earlystop/row_early_stop_process.ipp"
