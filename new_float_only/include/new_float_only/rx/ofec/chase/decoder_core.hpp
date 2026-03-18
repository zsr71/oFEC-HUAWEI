#pragma once

#include <vector>

#include "new_float_only/common/matrix/matrix.hpp"
#include "new_float_only/params.hpp"

namespace chase {

/**
 * 行级核心译码结果。
 * lout 存放每一行 256 维 extrinsic，produced_rows 标记该行是否真正产出结果。
 */
struct DecoderCoreResult {
  matrix::Matrix<float> lout;
  std::vector<bool> produced_rows;
};

/**
 * plain 行级译码核心。
 * 该函数逐行调用 hard decode 或 Chase decoder，并返回 256 维 extrinsic。
 */
DecoderCoreResult Decoder_Core_plain(const matrix::Matrix<float>& lin_matrix,
                                     const matrix::Matrix<float>& lch_matrix,
                                     bool use_hard_decode,
                                     const new_float_only::Params& config);

}  // namespace chase
