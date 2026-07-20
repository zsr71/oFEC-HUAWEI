#pragma once

namespace newcode {

template <typename LLR, typename Enable>
typename LinMatrixAdapter<LLR, Enable>::core_type
LinMatrixAdapter<LLR, Enable>::combine(const LLR& Lch, const LLR& La)
{
  // 输入:
  // - Lch: 信道 LLR。
  // - La : 先验/历史 LLR。
  // 输出:
  // - core_type 类型的合成结果。
  // 用途:
  // - 在“真实 LLR 域”中做线性合成，结果等价于 Lch + La，
  //   再按 core_type 的表示方式转回目标类型。
  const float sum = qfloat::llr_to_float(Lch) + qfloat::llr_to_float(La);
  return qfloat::llr_from_float<core_type>(sum);
}

template <typename LLR, typename Enable>
typename LinMatrixAdapter<LLR, Enable>::core_type
LinMatrixAdapter<LLR, Enable>::channel(const LLR& v)
{
  // 输入:
  // - v: 原始信道 LLR。
  // 输出:
  // - core_type 类型的信道值。
  // 用途:
  // - 仅保留信道项，不叠加先验信息；
  //   常用于单独构造 lch_matrix，供后续 Chase/Pyndiah 计算使用。
  return qfloat::llr_from_float<core_type>(qfloat::llr_to_float(v));
}

template <int NBITS, typename Store>
typename LinMatrixAdapter<qfloat::qfloat<NBITS, Store>>::core_type
LinMatrixAdapter<qfloat::qfloat<NBITS, Store>>::combine(const qfloat::qfloat<NBITS, Store>& Lch,
                                                const qfloat::qfloat<NBITS, Store>& La)
{
  // 输入:
  // - Lch: qfloat 形式的信道值。
  // - La : qfloat 形式的先验/历史值。
  // 输出:
  // - float 类型，但数值语义是“码值域”的和。
  // 用途:
  // - qfloat 特化路径不先反量化到真实幅度，而是直接把内部 code 相加，
  //   让 core 在量化码值域内工作，减少反复量化/反量化开销。
  return static_cast<float>(Lch.code() + La.code());
}

template <int NBITS, typename Store>
typename LinMatrixAdapter<qfloat::qfloat<NBITS, Store>>::core_type
LinMatrixAdapter<qfloat::qfloat<NBITS, Store>>::channel(const qfloat::qfloat<NBITS, Store>& v)
{
  // 输入:
  // - v: qfloat 形式的信道值。
  // 输出:
  // - float 类型的内部码值。
  // 用途:
  // - 为 qfloat 路径提供“仅信道项”的 core 表示，
  //   返回的是 code 域数值，不是反量化后的真实 LLR 幅度。
  return static_cast<float>(v.code());
}

template <typename LLR, typename Enable>
float ExtrinsicQuantizer<LLR, Enable>::quantize(float value)
{
  // 输入:
  // - value: 待写回的 extrinsic 浮点值。
  // 输出:
  // - 原样返回的 float。
  // 用途:
  // - 非量化类型路径不需要额外量化，直接保留浮点结果。
  return value;
}

template <int NBITS, typename Store>
float ExtrinsicQuantizer<qfloat::qfloat<NBITS, Store>>::quantize(float value)
{
  // 输入:
  // - value: 待写回的 extrinsic，当前按码值域 float 表示。
  // 输出:
  // - 经过 qfloat 饱和量化后，再转回 float 的结果。
  // 用途:
  // - 将 core 产生的 float 结果限制到 qfloat 可表达范围内，
  //   避免写回量化矩阵时出现越界。
  int code = static_cast<int>(std::lrint(value));
  const int lo = qfloat::qfloat<NBITS, Store>::LO();
  const int hi = qfloat::qfloat<NBITS, Store>::HI();
  if (code < lo) code = lo;
  if (code > hi) code = hi;
  qfloat::qfloat<NBITS, Store> q;
  q.set_code(code);
  return llr_to_float(q);
}

} // namespace newcode
