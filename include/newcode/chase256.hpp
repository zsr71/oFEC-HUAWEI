#pragma once
#include <cstddef>
#include <cstdint>
#include "newcode/params.hpp"

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

extern template void chase_decode_256_ebchPF<float >(const float*,  const float*,  float*,  const Params&);
extern template void chase_decode_256_ebchPF<int8_t>(const int8_t*, const int8_t*, int8_t*, const Params&);

} // namespace newcode
