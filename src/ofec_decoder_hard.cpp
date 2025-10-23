#include "newcode/ofec_decoder_hard.hpp"

#include "newcode/bch_255_239.hpp"
#include "newcode/llr_utils.hpp"

#include <array>
#include <cstdint>
#include <cmath>

namespace newcode {

template <typename LLR>
bool perform_hard_decode(const std::array<LLR, 256>& Lin256,
                         const std::array<LLR, 256>& Lch256,
                         std::array<LLR, 256>& Y2,
                         const Params& p)
{
  std::array<uint8_t, 256> hard_in{};
  for (int i = 0; i < 256; ++i) {
    hard_in[static_cast<size_t>(i)] = (llr_to_float(Lin256[static_cast<size_t>(i)]) < 0.f) ? 1u : 0u;
  }

  std::array<uint8_t, Params::BCH_N - 1> decoded{};
  if (!bch_255_239_decode_hiho_cw_255(hard_in.data(), decoded.data())) {
    return false;
  }

  std::array<uint8_t, Params::BCH_N> cw{};
  const int parity_len = static_cast<int>(Params::BCH_N) - 1;
  for (int i = 0; i < parity_len; ++i) {
    cw[static_cast<size_t>(i)] = decoded[static_cast<size_t>(i)];
  }

  uint8_t parity = 0;
  for (int i = 0; i < parity_len; ++i) {
    parity ^= cw[static_cast<size_t>(i)];
  }
  cw[static_cast<size_t>(Params::BCH_OVERALL_IDX)] = parity;

  const float hard_mag = std::fabs(p.HARD_LLR_MAG);
  for (int i = 0; i < static_cast<int>(Params::BCH_N); ++i) {
    const float sign = cw[static_cast<size_t>(i)] ? -1.f : 1.f;
    const float Lpost = sign * hard_mag;
    const float Lch = llr_to_float(Lch256[static_cast<size_t>(i)]);
    Y2[static_cast<size_t>(i)] = llr_from_float<LLR>(Lpost - Lch);
  }

  return true;
}

template bool perform_hard_decode<float>(const std::array<float, 256>&,
                                         const std::array<float, 256>&,
                                         std::array<float, 256>&,
                                         const Params&);
template bool perform_hard_decode<int8_t>(const std::array<int8_t, 256>&,
                                          const std::array<int8_t, 256>&,
                                          std::array<int8_t, 256>&,
                                          const Params&);
template bool perform_hard_decode<qfloat<4>>(const std::array<qfloat<4>, 256>&,
                                             const std::array<qfloat<4>, 256>&,
                                             std::array<qfloat<4>, 256>&,
                                             const Params&);
template bool perform_hard_decode<qfloat<5>>(const std::array<qfloat<5>, 256>&,
                                             const std::array<qfloat<5>, 256>&,
                                             std::array<qfloat<5>, 256>&,
                                             const Params&);

} // namespace newcode
