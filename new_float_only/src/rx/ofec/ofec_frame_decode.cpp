#include "new_float_only/ofec_decoder.hpp"
#include "new_float_only/params.hpp"
#include "new_float_only/rx/ofec/chase/decoder_core.hpp"

#include "detail/ofec_decode_impl.ipp"

namespace new_float_only {

/**
 * 顶层浮点 plain 解码入口。
 * 这里仅负责把公共参数传入 detail 层，真正的滑窗/Tile 逻辑在 .ipp 中实现。
 */
matrix::Matrix<float> decode_plain_llr(const matrix::Matrix<float>& channel_llr,
                                       const Params& config,
                                       bool normalize_extrinsic,
                                       const matrix::Matrix<float>* tx_llr_ref) {
  return detail::decode_plain_llr_impl(channel_llr,
                                       config,
                                       normalize_extrinsic,
                                       tx_llr_ref,
                                       &chase::Decoder_Core_plain);
}

}  // namespace new_float_only
