#include "newcode/decoder_api.hpp"

#include <iostream>
#include <mutex>

#include "newcode/llr_qpack.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {
namespace {

class EbchPFDecoder final : public IDecoder {
public:
  DecodeResult decode(const DecodeRequest& request) override {
    DecodeResult result;
    result.stats = {};

    switch (request.format) {
      case LlrFormat::Float:
        std::cout << "[INFO] (" << request.label << ") ebchPF decoder running in FLOAT\n";
        result.pre_decoder_llr = request.channel_llr;
        result.post_decoder_llr =
            ofec_decode_llr_ebchPF(result.pre_decoder_llr, request.params, &result.tile_stats,
                                   request.normalize_extrinsic);
        break;

      case LlrFormat::QFloat5: {
        std::cout << "[INFO] (" << request.label << ") ebchPF decoder running in qfloat<5>\n";
        auto quantized = quantize_matrix_to_qfloat<5>(request.channel_llr,
                                                      qfloat<5>::DEFAULT_CLIP);
        std::cout << "[INFO] (" << request.label << ") Quant clip=" << qfloat<5>::DEFAULT_CLIP
                  << " levels ±" << qfloat<5>::Q() << "\n";
        result.pre_decoder_llr = cast_matrix_from_qfloat(quantized);
        auto decoded = ofec_decode_llr_ebchPF(quantized, request.params, &result.tile_stats,
                                              request.normalize_extrinsic);
        result.post_decoder_llr = cast_matrix_from_qfloat(decoded);
        break;
      }

      case LlrFormat::QFloat4: {
        std::cout << "[INFO] (" << request.label << ") ebchPF decoder running in qfloat<4>\n";
        auto quantized = quantize_matrix_to_qfloat<4>(request.channel_llr,
                                                      qfloat<4>::DEFAULT_CLIP);
        std::cout << "[INFO] (" << request.label << ") Quant clip=" << qfloat<4>::DEFAULT_CLIP
                  << " levels ±" << qfloat<4>::Q() << "\n";
        result.pre_decoder_llr = cast_matrix_from_qfloat(quantized);
        auto decoded = ofec_decode_llr_ebchPF(quantized, request.params, &result.tile_stats,
                                              request.normalize_extrinsic);
        result.post_decoder_llr = cast_matrix_from_qfloat(decoded);
        break;
      }
    }

    result.stats.success = true;
    return result;
  }
};

static void ensure_registered() {
  static std::once_flag once;
  std::call_once(once, [] {
    register_decoder_factory([](const std::string& name) -> std::unique_ptr<IDecoder> {
      if (name == "ebchPF") {
        return std::make_unique<EbchPFDecoder>();
      }
      return nullptr;
    });
  });
}

} // namespace

void register_decoder_ebchPF_factory() {
  ensure_registered();
}

} // namespace newcode
