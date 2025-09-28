#pragma once
#include <cstddef>
#include <cstdint>
#include "newcode/params.hpp"

namespace newcode {

/**
 * @brief 单 Tile 内 256 维 Chase–Pyndiah：Y(256) -> extrinsic(256)
 *        仅流程 2~6（候选生成、BCH 硬译码、度量、竞争者、外信息）
 *
 * @tparam LLR  支持 float / int8_t / qfloat<NBITS>（要求：
 *              1) 可 static_cast<float>(LLR)；2) 若非算术类型，需提供 LLR::from_float(float)）
 *
 * @param Y256   输入 LLR（长度 256，第 255 号为整体奇偶位）
 * @param Y2_256 输出外信息 LLR（与输入同类型/同长度）
 * @param p      Params（使用 CHASE_* 与 CP_* 系数）
 */
template<typename LLR>
void chase_decode_256(const LLR* Y256, LLR* Y2_256, const Params& p);

// 显式实例化（与你项目中常用 LLR 类型对齐）
extern template void chase_decode_256<float >(const float*,  float*,  const Params&);
extern template void chase_decode_256<int8_t>(const int8_t*, int8_t*, const Params&);

} // namespace newcode
