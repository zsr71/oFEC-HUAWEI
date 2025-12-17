#include "newcode/ofec_decoder.hpp"
#include "newcode/decoder_core.hpp"

#include "detail/ofec_tile_impl.ipp"

namespace newcode {

template <typename LLR>
TileProcessResult<LLR> process_tile_plain(const Matrix<LLR>& tile_in,
                                          const Matrix<LLR>& ch_tile,
                                          const Params& p,
                                          size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic,
                                          const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return detail::process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                                   use_hard_decode, normalize_extrinsic,
                                   tx_llr_ref,
                                   &Decoder_Core_plain<CoreLLR>,
                                   /*last_tile_history_accum=*/nullptr,
                                   /*capture_last_tile_history=*/false);
}

template <typename LLR>
TileProcessResult<LLR> process_tile_ebchPF(const Matrix<LLR>& tile_in,
                                           const Matrix<LLR>& ch_tile,
                                           const Params& p,
                                           size_t tile_top_row_global,
                                           bool use_hard_decode,
                                           bool normalize_extrinsic,
                                           const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return detail::process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                                   use_hard_decode, normalize_extrinsic,
                                   tx_llr_ref,
                                   &Decoder_Core_ebchPF<CoreLLR>,
                                   /*last_tile_history_accum=*/nullptr,
                                   /*capture_last_tile_history=*/false);
}

// ===== 显式实例化 =====
template TileProcessResult<float>  process_tile_plain<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool, bool, const Matrix<float>*);
template TileProcessResult<float>  process_tile_ebchPF<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool, bool, const Matrix<float>*);

#define INSTANTIATE_TILE_QFLOAT(N) \
template TileProcessResult<newcode::qfloat<N>> process_tile_plain<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    const Params&, size_t, bool, bool, const Matrix<float>*); \
template TileProcessResult<newcode::qfloat<N>> process_tile_ebchPF<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    const Params&, size_t, bool, bool, const Matrix<float>*);

INSTANTIATE_TILE_QFLOAT(2)
INSTANTIATE_TILE_QFLOAT(3)
INSTANTIATE_TILE_QFLOAT(4)
INSTANTIATE_TILE_QFLOAT(5)
INSTANTIATE_TILE_QFLOAT(6)
INSTANTIATE_TILE_QFLOAT(7)
INSTANTIATE_TILE_QFLOAT(8)
INSTANTIATE_TILE_QFLOAT(9)
INSTANTIATE_TILE_QFLOAT(10)
INSTANTIATE_TILE_QFLOAT(11)
INSTANTIATE_TILE_QFLOAT(12)
INSTANTIATE_TILE_QFLOAT(13)
INSTANTIATE_TILE_QFLOAT(14)
INSTANTIATE_TILE_QFLOAT(15)

#undef INSTANTIATE_TILE_QFLOAT

} // namespace newcode
