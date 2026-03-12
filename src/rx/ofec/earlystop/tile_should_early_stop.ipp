#pragma once

namespace newcode {

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix)
{
  return tile_should_early_stop(lin_matrix, 1);
}

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix, int detect_mode)
{
  return detect_tile_early_stop(lin_matrix, detect_mode).all_rows_passed;
}

} // namespace newcode
