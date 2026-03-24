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
  // 输入:
  // - Lin256: 当前送入硬判 BCH 回退分支的 256 维输入 LLR。
  // - Lch256: 原始信道 LLR；当前实现未使用，但接口保留给统一流程。
  // - Y2    : 输出外信息/增量信息数组，由本函数填充。
  // - p     : 参数集合，主要使用 HARD_LLR_MAG 作为硬判后验幅度。
  // 输出:
  // - 返回 true 表示 BCH(255,239) 硬译码成功，Y2 已写入。
  // - 返回 false 表示硬译码失败，Y2 的内容不保证有效。
  // 用途:
  // - 当常规 Chase 软译码路径不可用或需要硬判回退时，
  //   用 Lin256 做硬判决 + BCH 硬译码，再构造一组简化的 extrinsic。
  (void)Lch256;

  std::array<uint8_t, 256> hard_in{};
  for (int i = 0; i < 256; ++i) {
    // 将输入 LLR 逐位硬判成 0/1。
    // 约定: LLR < 0 判为比特 1，否则判为比特 0。
    hard_in[static_cast<size_t>(i)] = (qfloat::llr_to_float(Lin256[static_cast<size_t>(i)]) < 0.f) ? 1u : 0u;
  }

  std::array<uint8_t, newcode::Params::BCH_N - 1> decoded{};
  if (!bch::bch_255_239_decode_hiho_cw_255(hard_in.data(), decoded.data())) {
    // BCH 硬译码失败则直接返回，让上层决定是否放弃该行/块。
    return false;
  }

  std::array<uint8_t, newcode::Params::BCH_N> cw{};
  const int parity_len = static_cast<int>(newcode::Params::BCH_N) - 1;
  for (int i = 0; i < parity_len; ++i) {
    // 将 BCH(255,239) 输出拷回到 256 长度码字前 255 位。
    cw[static_cast<size_t>(i)] = decoded[static_cast<size_t>(i)];
  }

  uint8_t parity = 0;
  for (int i = 0; i < parity_len; ++i) {
    // 重新计算整体奇偶校验位，补齐扩展 BCH 的第 256 位。
    parity ^= cw[static_cast<size_t>(i)];
  }
  cw[static_cast<size_t>(newcode::Params::BCH_OVERALL_IDX)] = parity;

  const float hard_mag = std::fabs(p.HARD_LLR_MAG);
  for (int i = 0; i < static_cast<int>(newcode::Params::BCH_N); ++i) {
    // 用一个固定幅度的硬判后验值表示译码结果:
    // - 译码比特为 0 -> +hard_mag
    // - 译码比特为 1 -> -hard_mag
    // 再减去当前输入 Lin256，形成写回上层的增量/外信息 Y2。
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
