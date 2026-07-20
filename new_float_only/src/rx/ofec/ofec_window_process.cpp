#include "new_float_only/ofec_decoder.hpp"
#include "new_float_only/params.hpp"
#include "new_float_only/rx/ofec/chase/decoder_core.hpp"

#include "detail/ofec_window_impl.ipp"

namespace new_float_only {

/**
 * 执行一个窗口内的全部 Tile。
 * 这个包装层主要负责把 public API 转交给 detail::process_window_impl。
 */
void process_window_plain(matrix::Matrix<float>& work_llr,
                          const matrix::Matrix<float>& channel_llr,
                          std::size_t win_start,
                          std::size_t win_end,
                          const Params& config,
                          std::size_t tile_height_rows,
                          std::size_t tile_stride_rows,
                          std::size_t tiles_per_window,
                          bool normalize_extrinsic,
                          const matrix::Matrix<float>* tx_llr_ref,
                          matrix::Matrix<float>* last_tile_history_accum) {
  detail::process_window_impl(work_llr,
                              channel_llr,
                              win_start,
                              win_end,
                              config,
                              tile_height_rows,
                              tile_stride_rows,
                              tiles_per_window,
                              normalize_extrinsic,
                              tx_llr_ref,
                              &chase::Decoder_Core_plain,
                              last_tile_history_accum);
}

}  // namespace new_float_only
