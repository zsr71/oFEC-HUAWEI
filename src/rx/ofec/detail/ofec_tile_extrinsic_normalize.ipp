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

