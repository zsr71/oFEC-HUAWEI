#pragma once

#include <cmath>

namespace newcode {

template <typename LLR>
TileEarlyStopResult detect_tile_early_stop_v2(const matrix::Matrix<LLR>& lin_matrix,
                                              const Params& p)
{
  TileEarlyStopResult res{};

  // 先初始化结果结构：记录总行数，并默认所有行都“不通过早停”。
  const size_t rows = lin_matrix.rows();
  const size_t cols = lin_matrix.cols();
  res.rows_total = rows;
  res.row_passed_flags.assign(rows, false);

  // 防御式检查：如果输入行宽不足 256，则认为当前 tile 形状异常，不触发早停。
  if (cols < 256) return res;

  // v2 判据的两个核心参数：
  // 1. |LLR| < threshold 视为“不可靠 bit”
  // 2. 一行中允许的不可靠 bit 数量上限为 max_unreliable_bits
  const float threshold = std::fabs(p.EARLY_STOP_V2_LLR_ABS_THRESHOLD);
  const int max_unreliable_bits = p.EARLY_STOP_V2_MAX_UNRELIABLE_BITS;

  bool all_pass = true;
  for (size_t r = 0; r < rows; ++r) {
    // 逐行统计“不可靠 bit”个数。
    int unreliable_count = 0;
    for (size_t j = 0; j < Params::BCH_N; ++j) {
      const float v = qfloat::llr_to_float(lin_matrix[r][j]);

      // 若该 bit 的 LLR 绝对值低于阈值，则认为它当前置信度不足。
      if (std::fabs(v) < threshold) {
        ++unreliable_count;

        // 一旦超出允许上限，就可以提前结束本行统计，
        // 因为这行已经不可能满足早停条件了。
        if (unreliable_count > max_unreliable_bits) {
          break;
        }
      }
    }

    // 本行不可靠 bit 数未超过阈值，则认为该行满足 v2 早停判据。
    const bool row_passed = unreliable_count <= max_unreliable_bits;

    if (row_passed) {
      ++res.rows_passed;
    } else {
      all_pass = false;
    }
    res.row_passed_flags[r] = row_passed;
  }

  // 只有当 tile 中每一行都满足 v2 判据时，才认为整个 tile 可以早停。
  res.all_rows_passed = all_pass && (rows > 0 ? res.rows_passed == rows : false);
  return res;
}

} // namespace newcode
