#include "newcode/ofec_decoder_hard.hpp"

#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"

#include <array>
#include <cstdint>
#include <cmath>

namespace newcode {

template <typename LLR>
bool perform_hard_decode(const std::array<LLR, 256>& Lin256,
                         const std::array<LLR, 256>& Lch256,
                         std::array<float, 256>& Y2,
                         const newcode::Params& p)
{
  std::array<uint8_t, 256> hard_in{};
  for (int i = 0; i < 256; ++i) {
    hard_in[static_cast<size_t>(i)] = (qfloat::llr_to_float(Lin256[static_cast<size_t>(i)]) < 0.f) ? 1u : 0u;
  }

  std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
  if (!bch::bch_255_239_decode_hiho_cw_255(hard_in.data(), decoded.data())) {
    return false;
  }

  std::array<uint8_t, newcode::Params::BCH_N> cw{};
  const int parity_len = static_cast<int>(newcode::Params::BCH_N) - 1;
  for (int i = 0; i < parity_len; ++i) {
    cw[static_cast<size_t>(i)] = decoded[static_cast<size_t>(i)];
  }

  uint8_t parity = 0;
  for (int i = 0; i < parity_len; ++i) {
    parity ^= cw[static_cast<size_t>(i)];
  }
  cw[static_cast<size_t>(newcode::Params::BCH_OVERALL_IDX)] = parity;

  const float hard_mag = std::fabs(p.HARD_LLR_MAG);
  for (int i = 0; i < static_cast<int>(newcode::Params::BCH_N); ++i) {
    const float sign = cw[static_cast<size_t>(i)] ? -1.f : 1.f;
    const float Lpost = sign * hard_mag;
    const float Lch = qfloat::llr_to_float(Lin256[static_cast<size_t>(i)]);
    Y2[static_cast<size_t>(i)] = Lpost - Lch;
  }

  return true;
}

template bool perform_hard_decode<float>(const std::array<float, 256>&,
                                         const std::array<float, 256>&,
                                         std::array<float, 256>&,
                                         const newcode::Params&);
template bool perform_hard_decode<int8_t>(const std::array<int8_t, 256>&,
                                          const std::array<int8_t, 256>&,
                                          std::array<float, 256>&,
                                          const newcode::Params&);

#define INSTANTIATE_HARD_DECODE_QFLOAT(N) \
template bool perform_hard_decode<qfloat::qfloat<N>>( \
    const std::array<qfloat::qfloat<N>, 256>&, \
    const std::array<qfloat::qfloat<N>, 256>&, \
    std::array<float, 256>&, \
    const newcode::Params&);

INSTANTIATE_HARD_DECODE_QFLOAT(2)
INSTANTIATE_HARD_DECODE_QFLOAT(3)
INSTANTIATE_HARD_DECODE_QFLOAT(4)
INSTANTIATE_HARD_DECODE_QFLOAT(5)
INSTANTIATE_HARD_DECODE_QFLOAT(6)
INSTANTIATE_HARD_DECODE_QFLOAT(7)
INSTANTIATE_HARD_DECODE_QFLOAT(8)
INSTANTIATE_HARD_DECODE_QFLOAT(9)
INSTANTIATE_HARD_DECODE_QFLOAT(10)
INSTANTIATE_HARD_DECODE_QFLOAT(11)
INSTANTIATE_HARD_DECODE_QFLOAT(12)
INSTANTIATE_HARD_DECODE_QFLOAT(13)
INSTANTIATE_HARD_DECODE_QFLOAT(14)
INSTANTIATE_HARD_DECODE_QFLOAT(15)

#undef INSTANTIATE_HARD_DECODE_QFLOAT

} // namespace newcode
