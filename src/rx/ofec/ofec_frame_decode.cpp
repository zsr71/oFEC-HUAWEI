#include "newcode/ofec_decoder.hpp"
#include "newcode/decoder_core.hpp"

#include "detail/ofec_decode_impl.ipp"

namespace newcode {

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_plain(const matrix::Matrix<LLR>& llr_mat, const Params& p,
                                  std::vector<TileEarlyStopCounter>* tile_stats,
                                  bool normalize_extrinsic,
                                  const matrix::Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return detail::ofec_decode_llr_impl(llr_mat, p, tile_stats, normalize_extrinsic,
                                      tx_llr_ref,
                                      &Decoder_Core_plain<CoreLLR>);
}

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_ebchPF(const matrix::Matrix<LLR>& llr_mat, const Params& p,
                                   std::vector<TileEarlyStopCounter>* tile_stats,
                                   bool normalize_extrinsic,
                                   const matrix::Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return detail::ofec_decode_llr_impl(llr_mat, p, tile_stats, normalize_extrinsic,
                                      tx_llr_ref,
                                      &Decoder_Core_ebchPF<CoreLLR>);
}

// ===== 显式实例化 =====
template matrix::Matrix<float>  ofec_decode_llr_plain<float >(const matrix::Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);
template matrix::Matrix<float>  ofec_decode_llr_ebchPF<float >(const matrix::Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);

#define INSTANTIATE_DECODE_QFLOAT(N) \
template matrix::Matrix<qfloat::qfloat<N>> ofec_decode_llr_plain<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>&, const Params&, \
    std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*); \
template matrix::Matrix<qfloat::qfloat<N>> ofec_decode_llr_ebchPF<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>&, const Params&, \
    std::vector<TileEarlyStopCounter>*, bool, const matrix::Matrix<float>*);
INSTANTIATE_DECODE_QFLOAT(2)
INSTANTIATE_DECODE_QFLOAT(3)
INSTANTIATE_DECODE_QFLOAT(4)
INSTANTIATE_DECODE_QFLOAT(5)
INSTANTIATE_DECODE_QFLOAT(6)
INSTANTIATE_DECODE_QFLOAT(7)
INSTANTIATE_DECODE_QFLOAT(8)
INSTANTIATE_DECODE_QFLOAT(9)
INSTANTIATE_DECODE_QFLOAT(10)
INSTANTIATE_DECODE_QFLOAT(11)
INSTANTIATE_DECODE_QFLOAT(12)
INSTANTIATE_DECODE_QFLOAT(13)
INSTANTIATE_DECODE_QFLOAT(14)
INSTANTIATE_DECODE_QFLOAT(15)

#undef INSTANTIATE_DECODE_QFLOAT

} // namespace newcode
