#pragma once

namespace newcode {

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix)
{
  return tile_should_early_stop(lin_matrix, Params{}, 1);
}

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix, int detect_mode)
{
  return tile_should_early_stop(lin_matrix, Params{}, detect_mode);
}

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix,
                            const Params& p,
                            int detect_mode)
{
  return detect_tile_early_stop(lin_matrix, p, detect_mode).all_rows_passed;
}

} // namespace newcode
