#include "new_float_only/ofec_decoder.hpp"
#include "new_float_only/params.hpp"
#include "new_float_only/rx/ofec/chase/decoder_core.hpp"

#include "detail/ofec_tile_impl.ipp"

namespace new_float_only {

/**
 * 单个 Tile 的 public 包装层。
 * 这里不再做 MUX/early-stop 选择，所有行都会直接进入 plain core。
 */
TileProcessResult process_tile_plain(const matrix::Matrix<float>& tile_in,
                                     const matrix::Matrix<float>& ch_tile,
                                     const Params& config,
                                     std::size_t tile_top_row_global,
                                     bool use_hard_decode,
                                     bool normalize_extrinsic,
                                     const matrix::Matrix<float>* tx_llr_ref,
                                     matrix::Matrix<float>* last_tile_history_accum,
                                     bool capture_last_tile_history) {
  return detail::process_tile_impl(tile_in,
                                   ch_tile,
                                   config,
                                   tile_top_row_global,
                                   use_hard_decode,
                                   normalize_extrinsic,
                                   tx_llr_ref,
                                   &chase::Decoder_Core_plain,
                                   last_tile_history_accum,
                                   capture_last_tile_history);
}

}  // namespace new_float_only
