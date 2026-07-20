#pragma once

namespace newcode {

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix)
{
  return detect_tile_early_stop(lin_matrix, Params{}).all_rows_passed;
}

template <typename LLR>
bool tile_should_early_stop(const Matrix<LLR>& lin_matrix,
                            const Params& p)
{
  return detect_tile_early_stop(lin_matrix, p).all_rows_passed;
}

} // namespace newcode
