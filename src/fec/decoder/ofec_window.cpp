#include "newcode/ofec_decoder.hpp"
#include "newcode/decoder_core.hpp"

#include "detail/ofec_window_impl.ipp"

namespace newcode {

template <typename LLR>
void process_window_plain(Matrix<LLR>& work_llr,
                          const Matrix<LLR>& channel_llr,
                          size_t win_start, size_t win_end, const Params& p,
                          size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                          std::vector<TileEarlyStopCounter>* tile_stats,
                          bool normalize_extrinsic,
                          const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  detail::process_window_impl(work_llr, channel_llr, win_start, win_end, p,
                              tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                              tile_stats, normalize_extrinsic, tx_llr_ref,
                              &Decoder_Core_plain<CoreLLR>,
                              /*last_tile_history_accum=*/nullptr);
}

template <typename LLR>
void process_window_ebchPF(Matrix<LLR>& work_llr,
                           const Matrix<LLR>& channel_llr,
                           size_t win_start, size_t win_end, const Params& p,
                           size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                           std::vector<TileEarlyStopCounter>* tile_stats,
                           bool normalize_extrinsic,
                           const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  detail::process_window_impl(work_llr, channel_llr, win_start, win_end, p,
                              tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                              tile_stats, normalize_extrinsic, tx_llr_ref,
                              &Decoder_Core_ebchPF<CoreLLR>,
                              /*last_tile_history_accum=*/nullptr);
}

// ===== 显式实例化 =====
template void process_window_plain<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool,
                                           const Matrix<float>*);

template void process_window_ebchPF<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool,
                                           const Matrix<float>*);


template void process_window<float >(Matrix<float>&,  const Matrix<float>&,
                                     std::size_t, std::size_t, const Params&,
                                     std::size_t, std::size_t, std::size_t,
                                     std::vector<TileEarlyStopCounter>*, bool,
                                     const Matrix<float>*);


#define INSTANTIATE_WINDOW_QFLOAT(N) \
template void process_window_plain<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    size_t, size_t, const Params&, \
    size_t, size_t, size_t, \
    std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*); \
template void process_window_ebchPF<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    size_t, size_t, const Params&, \
    size_t, size_t, size_t, \
    std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*); \
template void process_window<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    std::size_t, std::size_t, const Params&, \
    std::size_t, std::size_t, std::size_t, \
    std::vector<TileEarlyStopCounter>*, bool, \
    const Matrix<float>*);

INSTANTIATE_WINDOW_QFLOAT(2)
INSTANTIATE_WINDOW_QFLOAT(3)
INSTANTIATE_WINDOW_QFLOAT(4)
INSTANTIATE_WINDOW_QFLOAT(5)
INSTANTIATE_WINDOW_QFLOAT(6)
INSTANTIATE_WINDOW_QFLOAT(7)
INSTANTIATE_WINDOW_QFLOAT(8)
INSTANTIATE_WINDOW_QFLOAT(9)
INSTANTIATE_WINDOW_QFLOAT(10)
INSTANTIATE_WINDOW_QFLOAT(11)
INSTANTIATE_WINDOW_QFLOAT(12)
INSTANTIATE_WINDOW_QFLOAT(13)
INSTANTIATE_WINDOW_QFLOAT(14)
INSTANTIATE_WINDOW_QFLOAT(15)

#undef INSTANTIATE_WINDOW_QFLOAT

} // namespace newcode
