#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace newcode {
namespace detail {

template <typename Float>
void normalize_extrinsic_lout(matrix::Matrix<Float>& lout,
                              const std::vector<bool>& produced_rows,
                              Float beta)
{
  // 输入:
  // - lout: decoder core 输出的外信息矩阵，会被原地缩放。
  // - produced_rows: 哪些 decoder row 实际有输出。
  // - beta: Pyndiah fallback 幅度，用于识别“退化外信息”。
  // 输出:
  // - 无返回值；直接修改 lout。
  // 用途:
  // - 估计当前 tile 外信息的平均幅度，并统一缩放到更稳定的数值范围，
  //   减少不同 tile/轮次之间幅度漂移。
  auto is_fallback = [&](Float w) -> bool {
    // 当 |w| 接近 beta 时，说明它可能来自 fallback 路径；
    // 这类值用 |w/beta| 参与统计，避免把固定幅度直接当作真实置信度。
    const Float target = beta;
    const Float diff   = std::fabs(std::fabs(w) - target);
    const Float tol    = static_cast<Float>(1e-4f) * std::max(static_cast<Float>(1.0f), target);
    return diff <= tol;
  };

  double acc = 0.0;
  std::size_t cnt = 0;

  const std::size_t Rcnt = lout.rows();
  const std::size_t Ccnt = lout.cols();

  for (std::size_t r = 0; r < Rcnt; ++r)
  {
    if (r >= produced_rows.size() || !produced_rows[r]) continue;
    for (std::size_t j = 0; j < Ccnt; ++j)
    {
      const Float w = lout[r][j];
      if (is_fallback(w)) {acc += std::fabs(w/beta); ++cnt; continue;};
      acc += std::fabs(w);
      ++cnt;
    }
  }

  if (cnt > 0)
  {
    const Float g_alpha = static_cast<Float>(acc / static_cast<double>(cnt));
    if (g_alpha > static_cast<Float>(0.0f))
    {
      // 用平均绝对值的倒数作为归一化系数，把当前 tile 的外信息整体拉回统一尺度。
      const Float scale = static_cast<Float>(1.0f) / g_alpha;
      for (std::size_t r = 0; r < Rcnt; ++r)
      {
        if (r >= produced_rows.size() || !produced_rows[r]) continue;
        for (std::size_t j = 0; j < Ccnt; ++j)
        {
          lout[r][j] *= scale;
        }
      }
    }
  }
}

} // namespace detail
} // namespace newcode
