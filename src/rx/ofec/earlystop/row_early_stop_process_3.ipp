#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/params.hpp"

namespace newcode {

template <typename LLR>
bool row_early_stop_process_3(const LLR* lin256,
                              const LLR* lch256,
                              float* y2_256,
                              const newcode::Params& p)
{
  (void)lch256;

  std::array<uint8_t, newcode::Params::BCH_N> hard_in{};
  for (size_t j = 0; j < newcode::Params::BCH_N; ++j) {
    // mode3 的硬判完全基于当前 lin256，不再叠加额外软信息修正。
    hard_in[j] = (qfloat::llr_to_float(lin256[j]) < 0.0f) ? 1u : 0u;
  }

  std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
  if (!bch::bch_255_239_decode_hiho_cw_255(hard_in.data(), decoded.data())) {
    return false;
  }

  std::array<uint8_t, newcode::Params::BCH_N> cw{};
  uint8_t overall = 0u;
  for (size_t j = 0; j < newcode::Params::BCH_N - 1; ++j) {
    cw[j] = decoded[j];
    overall ^= cw[j];
  }
  cw[static_cast<size_t>(newcode::Params::BCH_OVERALL_IDX)] = overall;

  const float hard_mag = std::fabs(p.EARLY_STOP_ACTION_HARD_LLR_MAG);
  for (size_t j = 0; j < newcode::Params::BCH_N; ++j) {
    y2_256[j] = cw[j] ? -hard_mag : hard_mag;
  }
  return true;
}

}  // namespace newcode
