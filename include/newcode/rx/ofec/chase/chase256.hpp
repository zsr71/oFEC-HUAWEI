#pragma once
#include <cstddef>
#include <cstdint>
#include "newcode/params.hpp"
#include "newcode/common/qfloat/qfloat.hpp"

namespace chase {

template<typename LLR>
void chase_decode_256_plain(const LLR* Lin256,
                            const LLR* Lch256,
                            float* Y2_256,
                            const newcode::Params& p);

template<typename LLR>
void chase_decode_256_plain(const LLR* Y256,
                            float* Y2_256,
                            const newcode::Params& p);

template<typename LLR>
void chase_decode_256_ebchPF(const LLR* Lin256,
                             const LLR* Lch256,
                             float* Y2_256,
                             const newcode::Params& p);

template<typename LLR>
void chase_decode_256_ebchPF(const LLR* Y256,
                             float* Y2_256,
                             const newcode::Params& p);

template<typename LLR>
void chase_decode_256_topk_pruned(const LLR* Lin256,
                                  const LLR* Lch256,
                                  float* Y2_256,
                                  const newcode::Params& p);

template<typename LLR>
void chase_decode_256_topk_pruned(const LLR* Y256,
                                  float* Y2_256,
                                  const newcode::Params& p);

template<typename LLR>
void chase_decode_256_global_pair(const LLR* Lin256,
                                  const LLR* Lch256,
                                  float* Y2_256,
                                  const newcode::Params& p);

template<typename LLR>
void chase_decode_256_global_pair(const LLR* Y256,
                                  float* Y2_256,
                                  const newcode::Params& p);

template<typename LLR>
void chase_decode_256_group_minima(const LLR* Lin256,
                                   const LLR* Lch256,
                                   float* Y2_256,
                                   const newcode::Params& p);

template<typename LLR>
void chase_decode_256_group_minima(const LLR* Y256,
                                   float* Y2_256,
                                   const newcode::Params& p);

// 显式实例化（与你项目中常用 LLR 类型对齐）
extern template void chase_decode_256_plain<float >(const float*,  const float*,  float*,  const newcode::Params&);
extern template void chase_decode_256_plain<float >(const float*,  float*,  const newcode::Params&);    
extern template void chase_decode_256_topk_pruned<float >(const float*,  const float*,  float*,  const newcode::Params&);
extern template void chase_decode_256_topk_pruned<float >(const float*,  float*,  const newcode::Params&);
extern template void chase_decode_256_global_pair<float >(const float*,  const float*,  float*,  const newcode::Params&);
extern template void chase_decode_256_global_pair<float >(const float*,  float*,  const newcode::Params&);
extern template void chase_decode_256_group_minima<float >(const float*,  const float*,  float*,  const newcode::Params&);
extern template void chase_decode_256_group_minima<float >(const float*,  float*,  const newcode::Params&);
#define DECLARE_CHASE256_QFLOAT(N) \
extern template void chase_decode_256_plain<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_plain<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_ebchPF<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_ebchPF<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_topk_pruned<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_topk_pruned<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_global_pair<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_global_pair<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_group_minima<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, const qfloat::qfloat<N>*, float*, const newcode::Params&); \
extern template void chase_decode_256_group_minima<qfloat::qfloat<N>>(const qfloat::qfloat<N>*, float*, const newcode::Params&);

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

extern template void chase_decode_256_ebchPF<float >(const float*,  const float*,  float*,  const newcode::Params&);
extern template void chase_decode_256_ebchPF<float >(const float*,  float*,  const newcode::Params&);

} // namespace newcode
