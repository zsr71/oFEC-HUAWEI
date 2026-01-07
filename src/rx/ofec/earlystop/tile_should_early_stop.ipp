#pragma once

namespace newcode {

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix)
{

  return tile_early_stop_stats1(lin_matrix).all_rows_passed;
}

} // namespace newcode
