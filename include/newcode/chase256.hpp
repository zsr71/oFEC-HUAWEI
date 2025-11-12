#pragma once
#include <cstddef>
#include <cstdint>
#include "newcode/params.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {

template<typename LLR>
void chase_decode_256_plain(const LLR* Lin256,
                            const LLR* Lch256,
                            LLR* Y2_256,
                            const Params& p);

template<typename LLR>
void chase_decode_256_plain(const LLR* Y256,
                            LLR* Y2_256,
                            const Params& p);

template<typename LLR>
void chase_decode_256_ebchPF(const LLR* Lin256,
                             const LLR* Lch256,
                             LLR* Y2_256,
                             const Params& p);

template<typename LLR>
void chase_decode_256_ebchPF(const LLR* Y256,
                             LLR* Y2_256,
                             const Params& p);

// 显式实例化（与你项目中常用 LLR 类型对齐）
extern template void chase_decode_256_plain<float >(const float*,  const float*,  float*,  const Params&);
extern template void chase_decode_256_plain<int8_t>(const int8_t*, const int8_t*, int8_t*, const Params&);
extern template void chase_decode_256_plain<float >(const float*,  float*,  const Params&);
extern template void chase_decode_256_plain<int8_t>(const int8_t*, int8_t*, const Params&);

#define DECLARE_CHASE256_QFLOAT(N) \
extern template void chase_decode_256_plain<qfloat<N>>(const qfloat<N>*, const qfloat<N>*, qfloat<N>*, const Params&); \
extern template void chase_decode_256_plain<qfloat<N>>(const qfloat<N>*, qfloat<N>*, const Params&); \
extern template void chase_decode_256_ebchPF<qfloat<N>>(const qfloat<N>*, const qfloat<N>*, qfloat<N>*, const Params&); \
extern template void chase_decode_256_ebchPF<qfloat<N>>(const qfloat<N>*, qfloat<N>*, const Params&);

DECLARE_CHASE256_QFLOAT(2)
DECLARE_CHASE256_QFLOAT(3)
DECLARE_CHASE256_QFLOAT(4)
DECLARE_CHASE256_QFLOAT(5)
DECLARE_CHASE256_QFLOAT(6)
DECLARE_CHASE256_QFLOAT(7)
DECLARE_CHASE256_QFLOAT(8)
DECLARE_CHASE256_QFLOAT(9)
DECLARE_CHASE256_QFLOAT(10)
DECLARE_CHASE256_QFLOAT(11)
DECLARE_CHASE256_QFLOAT(12)
DECLARE_CHASE256_QFLOAT(13)
DECLARE_CHASE256_QFLOAT(14)
DECLARE_CHASE256_QFLOAT(15)

#undef DECLARE_CHASE256_QFLOAT

extern template void chase_decode_256_ebchPF<float >(const float*,  const float*,  float*,  const Params&);
extern template void chase_decode_256_ebchPF<int8_t>(const int8_t*, const int8_t*, int8_t*, const Params&);
extern template void chase_decode_256_ebchPF<float >(const float*,  float*,  const Params&);
extern template void chase_decode_256_ebchPF<int8_t>(const int8_t*, int8_t*, const Params&);

} // namespace newcode
