#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace new_float_only {
namespace detail {

// 对当前 tile 产生的外信息矩阵做整体幅度归一化。
//
// 这里的目标不是改变符号方向，而是把本轮真正产出的 extrinsic 调整到一个更稳定
// 的平均绝对值尺度，避免不同 tile / 不同参数下外信息幅度漂移过大。
//
// 归一化时只统计 produced_rows 中为真的那些行；未产出输出的行保持原样。
template <typename Float>
void normalize_extrinsic_lout(matrix::Matrix<Float>& lout,
                              const std::vector<bool>& produced_rows,
                              Float beta)
{
  // 某些位置的 omega 可能直接走了 fallback，幅度会恰好落在 ±beta。
  // 这些值在统计平均绝对值时单独按 |w / beta| 计入，避免 beta 选得很大时把
  // 均值尺度直接拉爆。
  auto is_fallback = [&](Float w) -> bool {
    const Float target = beta;
    const Float diff   = std::fabs(std::fabs(w) - target);
    const Float tol    = static_cast<Float>(1e-4f) * std::max(static_cast<Float>(1.0f), target);
    return diff <= tol;
  };

  double acc = 0.0;
  std::size_t cnt = 0;

  const std::size_t Rcnt = lout.rows();
  const std::size_t Ccnt = lout.cols();

  // 先统计本轮有效输出的平均绝对值。只有 produced 的行才参与统计。
  for (std::size_t r = 0; r < Rcnt; ++r)
  {
    if (r >= produced_rows.size() || !produced_rows[r]) continue;
    for (std::size_t j = 0; j < Ccnt; ++j)
    {
      const Float w = lout[r][j];
      // 如果这个位置是 fallback 直接给出的 ±beta，就按 |w / beta| 计入。
      // 这样处理后，fallback 位对均值的贡献大致归一到 1，而不是随着 beta 的
      // 绝对值线性放大，避免 beta 本身主导整轮归一化尺度。
      if (is_fallback(w)) {continue;};
      // 非 fallback 的正常 extrinsic 直接按自身绝对值参与均值统计。
      acc += std::fabs(w);
      ++cnt;
    }
  }

  if (cnt > 0)
  {
    // g_alpha 表示当前 tile 外信息的平均绝对值；后面用 1 / g_alpha 做统一缩放，
    // 让本轮输出整体回到“平均绝对值约为 1”的量级。
    const Float g_alpha = static_cast<Float>(acc / static_cast<double>(cnt));
    if (g_alpha > static_cast<Float>(0.0f))
    {
      const Float scale = static_cast<Float>(1.0f) / g_alpha;
      for (std::size_t r = 0; r < Rcnt; ++r)
      {
        if (r >= produced_rows.size() || !produced_rows[r]) continue;
        for (std::size_t j = 0; j < Ccnt; ++j)
        {
          // 归一化只改幅度，不改符号。
          lout[r][j] *= scale;
        }
      }
    }
  }
}

} // namespace detail
} // namespace new_float_only
