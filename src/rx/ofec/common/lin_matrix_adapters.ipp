#pragma once

namespace newcode {

template <typename LLR, typename Enable>
typename LinMatrixAdapter<LLR, Enable>::core_type
LinMatrixAdapter<LLR, Enable>::combine(const LLR& Lch, const LLR& La)
{
  const float sum = llr_to_float(Lch) + llr_to_float(La);
  return llr_from_float<core_type>(sum);
}

template <typename LLR, typename Enable>
typename LinMatrixAdapter<LLR, Enable>::core_type
LinMatrixAdapter<LLR, Enable>::channel(const LLR& v)
{
  return llr_from_float<core_type>(llr_to_float(v));
}

template <int NBITS, typename Store>
typename LinMatrixAdapter<qfloat<NBITS, Store>>::core_type
LinMatrixAdapter<qfloat<NBITS, Store>>::combine(const qfloat<NBITS, Store>& Lch,
                                                const qfloat<NBITS, Store>& La)
{
  return static_cast<float>(Lch.code() + La.code());
}

template <int NBITS, typename Store>
typename LinMatrixAdapter<qfloat<NBITS, Store>>::core_type
LinMatrixAdapter<qfloat<NBITS, Store>>::channel(const qfloat<NBITS, Store>& v)
{
  return static_cast<float>(v.code());
}

template <typename LLR, typename Enable>
float ExtrinsicQuantizer<LLR, Enable>::quantize(float value)
{
  return value;
}

template <int NBITS, typename Store>
float ExtrinsicQuantizer<qfloat<NBITS, Store>>::quantize(float value)
{
  int code = static_cast<int>(std::lrint(value));
  const int lo = qfloat<NBITS, Store>::LO();
  const int hi = qfloat<NBITS, Store>::HI();
  if (code < lo) code = lo;
  if (code > hi) code = hi;
  qfloat<NBITS, Store> q;
  q.set_code(code);
  return llr_to_float(q);
}

} // namespace newcode
