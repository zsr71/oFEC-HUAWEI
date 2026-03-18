#include "new_float_only/ofec_decoder_hard.hpp"

#include "new_float_only/common/bch/bch_255_239.hpp"

#include <array>
#include <cmath>

namespace new_float_only {

/**
 * 硬判决译码。
 * 这里刻意保留旧 plain 流程的输出口径：最终 extrinsic 仍减去 lin256，
 * 而不是减去 lch256，从而保持与旧库当前行为一致。
 */
bool perform_hard_decode(const std::array<float, 256>& lin256,
                         const std::array<float, 256>& lch256,
                         std::array<float, 256>& y2,
                         const new_float_only::Params& p) {
  (void)lch256;
  std::array<uint8_t, 256> hard_in{};
  for (int i = 0; i < 256; ++i) {
    hard_in[static_cast<std::size_t>(i)] = (lin256[static_cast<std::size_t>(i)] < 0.f) ? 1u : 0u;
  }

  std::array<uint8_t, new_float_only::Params::BCH_N - 1> decoded{};
  if (!bch::bch_255_239_decode_hiho_cw_255(hard_in.data(), decoded.data())) {
    return false;
  }

  std::array<uint8_t, new_float_only::Params::BCH_N> cw{};
  const int parity_len = static_cast<int>(new_float_only::Params::BCH_N) - 1;
  for (int i = 0; i < parity_len; ++i) {
    cw[static_cast<size_t>(i)] = decoded[static_cast<size_t>(i)];
  }

  uint8_t parity = 0;
  for (int i = 0; i < parity_len; ++i) {
    parity ^= cw[static_cast<size_t>(i)];
  }
  cw[static_cast<size_t>(new_float_only::Params::BCH_OVERALL_IDX)] = parity;

  const float hard_mag = std::fabs(p.HARD_LLR_MAG);
  for (int i = 0; i < static_cast<int>(new_float_only::Params::BCH_N); ++i) {
    const float sign = cw[static_cast<std::size_t>(i)] ? -1.f : 1.f;
    const float Lpost = sign * hard_mag;
    const float Lin = lin256[static_cast<std::size_t>(i)];
    y2[static_cast<std::size_t>(i)] = Lpost - Lin;
  }

  return true;
}

}  // namespace new_float_only
