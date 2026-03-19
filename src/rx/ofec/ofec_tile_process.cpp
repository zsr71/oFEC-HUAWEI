#include "newcode/ofec_decoder.hpp"
#include "newcode/params.hpp"
#include "newcode/ofec/mux/mux_siso_budget.hpp"
#include "newcode/rx/ofec/chase/decoder_core.hpp"

#include "detail/ofec_tile_impl.ipp"

namespace newcode {

template <typename LLR>
TileProcessResult<LLR> process_tile_plain(const matrix::Matrix<LLR>& tile_in,
                                          const matrix::Matrix<LLR>& ch_tile,
                                          const newcode::Params& p,
                                          size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic,
                                          const matrix::Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  const int siso_active_for_tile =
      newcode::mux::pick_siso_active_for_tile(p.SISO_ACTIVE_LIST, 0);
  return detail::process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                                   siso_active_for_tile,
                                   use_hard_decode, normalize_extrinsic,
                                   tx_llr_ref,
                                   &chase::Decoder_Core_plain<CoreLLR>,
                                   /*last_tile_history_accum=*/nullptr,
                                   /*capture_last_tile_history=*/false);
}

template <typename LLR>
TileProcessResult<LLR> process_tile_ebchPF(const matrix::Matrix<LLR>& tile_in,
                                           const matrix::Matrix<LLR>& ch_tile,
                                           const newcode::Params& p,
                                           size_t tile_top_row_global,
                                           bool use_hard_decode,
                                           bool normalize_extrinsic,
                                           const matrix::Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  const int siso_active_for_tile =
      newcode::mux::pick_siso_active_for_tile(p.SISO_ACTIVE_LIST, 0);
  return detail::process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                                   siso_active_for_tile,
                                   use_hard_decode, normalize_extrinsic,
                                   tx_llr_ref,
                                   &chase::Decoder_Core_ebchPF<CoreLLR>,
                                   /*last_tile_history_accum=*/nullptr,
                                   /*capture_last_tile_history=*/false);
}

template <typename LLR>
TileProcessResult<LLR> process_tile_topk_pruned(const matrix::Matrix<LLR>& tile_in,
                                                const matrix::Matrix<LLR>& ch_tile,
                                                const newcode::Params& p,
                                                size_t tile_top_row_global,
                                                bool use_hard_decode,
                                                bool normalize_extrinsic,
                                                const matrix::Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  const int siso_active_for_tile =
      newcode::mux::pick_siso_active_for_tile(p.SISO_ACTIVE_LIST, 0);
  return detail::process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                                   siso_active_for_tile,
                                   use_hard_decode, normalize_extrinsic,
                                   tx_llr_ref,
                                   &chase::Decoder_Core_topk_pruned<CoreLLR>,
                                   /*last_tile_history_accum=*/nullptr,
                                   /*capture_last_tile_history=*/false);
}

template <typename LLR>
TileProcessResult<LLR> process_tile_global_pair(const matrix::Matrix<LLR>& tile_in,
                                                const matrix::Matrix<LLR>& ch_tile,
                                                const newcode::Params& p,
                                                size_t tile_top_row_global,
                                                bool use_hard_decode,
                                                bool normalize_extrinsic,
                                                const matrix::Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  const int siso_active_for_tile =
      newcode::mux::pick_siso_active_for_tile(p.SISO_ACTIVE_LIST, 0);
  return detail::process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                                   siso_active_for_tile,
                                   use_hard_decode, normalize_extrinsic,
                                   tx_llr_ref,
                                   &chase::Decoder_Core_global_pair<CoreLLR>,
                                   /*last_tile_history_accum=*/nullptr,
                                   /*capture_last_tile_history=*/false);
}

template <typename LLR>
TileProcessResult<LLR> process_tile_group_minima(const matrix::Matrix<LLR>& tile_in,
                                                 const matrix::Matrix<LLR>& ch_tile,
                                                 const newcode::Params& p,
                                                 size_t tile_top_row_global,
                                                 bool use_hard_decode,
                                                 bool normalize_extrinsic,
                                                 const matrix::Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  const int siso_active_for_tile =
      newcode::mux::pick_siso_active_for_tile(p.SISO_ACTIVE_LIST, 0);
  return detail::process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                                   siso_active_for_tile,
                                   use_hard_decode, normalize_extrinsic,
                                   tx_llr_ref,
                                   &chase::Decoder_Core_group_minima<CoreLLR>,
                                   /*last_tile_history_accum=*/nullptr,
                                   /*capture_last_tile_history=*/false);
}

// ===== 显式实例化 =====
template TileProcessResult<float>  process_tile_plain<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,  const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*);
template TileProcessResult<float>  process_tile_ebchPF<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,  const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*);
template TileProcessResult<float>  process_tile_topk_pruned<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,  const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*);
template TileProcessResult<float>  process_tile_global_pair<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,  const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*);
template TileProcessResult<float>  process_tile_group_minima<float >(const matrix::Matrix<float>&,  const matrix::Matrix<float>&,  const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*);

#define INSTANTIATE_TILE_QFLOAT(N) \
template TileProcessResult<qfloat::qfloat<N>> process_tile_plain<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>&, const matrix::Matrix<qfloat::qfloat<N>>&, \
    const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*); \
template TileProcessResult<qfloat::qfloat<N>> process_tile_ebchPF<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>&, const matrix::Matrix<qfloat::qfloat<N>>&, \
    const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*); \
template TileProcessResult<qfloat::qfloat<N>> process_tile_topk_pruned<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>&, const matrix::Matrix<qfloat::qfloat<N>>&, \
    const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*); \
template TileProcessResult<qfloat::qfloat<N>> process_tile_global_pair<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>&, const matrix::Matrix<qfloat::qfloat<N>>&, \
    const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*); \
template TileProcessResult<qfloat::qfloat<N>> process_tile_group_minima<qfloat::qfloat<N>>( \
    const matrix::Matrix<qfloat::qfloat<N>>&, const matrix::Matrix<qfloat::qfloat<N>>&, \
    const newcode::Params&, size_t, bool, bool, const matrix::Matrix<float>*);
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
