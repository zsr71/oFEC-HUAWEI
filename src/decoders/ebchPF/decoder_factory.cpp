#include "newcode/decoder_api.hpp"

#include <iostream>
#include <mutex>
#include <stdexcept>

#include "newcode/llr_qpack.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {
namespace {

template <int NBITS>
void decode_ebchpf_qfloat(const DecodeRequest& request, DecodeResult& result) {
  using Q = qfloat<NBITS>;
  const float clip = (request.quant_clip > 0.0f) ? request.quant_clip : Q::DEFAULT_CLIP;
  std::cout << "[INFO] (" << request.label << ") ebchPF decoder running in qfloat<"
            << NBITS << ">\n";
  auto quantized = quantize_matrix_to_qfloat<NBITS>(request.channel_llr, clip);
  std::cout << "[INFO] (" << request.label << ") Quant clip=" << clip
            << " levels ±" << Q::Q() << "\n";
  result.pre_decoder_llr = dequantize_matrix_from_qfloat(quantized, clip);
  auto decoded = ofec_decode_llr_ebchPF(quantized, request.params, &result.tile_stats,
                                        request.normalize_extrinsic);
  result.post_decoder_llr = dequantize_matrix_from_qfloat(decoded, clip);
}

void decode_ebchpf_quantized(const DecodeRequest& request, DecodeResult& result) {
  switch (request.quant_bits) {
    case 2:  decode_ebchpf_qfloat<2 >(request, result); break;
    case 3:  decode_ebchpf_qfloat<3 >(request, result); break;
    case 4:  decode_ebchpf_qfloat<4 >(request, result); break;
    case 5:  decode_ebchpf_qfloat<5 >(request, result); break;
    case 6:  decode_ebchpf_qfloat<6 >(request, result); break;
    case 7:  decode_ebchpf_qfloat<7 >(request, result); break;
    case 8:  decode_ebchpf_qfloat<8 >(request, result); break;
    case 9:  decode_ebchpf_qfloat<9 >(request, result); break;
    case 10: decode_ebchpf_qfloat<10>(request, result); break;
    case 11: decode_ebchpf_qfloat<11>(request, result); break;
    case 12: decode_ebchpf_qfloat<12>(request, result); break;
    case 13: decode_ebchpf_qfloat<13>(request, result); break;
    case 14: decode_ebchpf_qfloat<14>(request, result); break;
    case 15: decode_ebchpf_qfloat<15>(request, result); break;
    default:
      throw std::runtime_error("[ERROR] Unsupported quant_bits for ebchPF decoder: " +
                               std::to_string(request.quant_bits));
  }
}

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

      case LlrFormat::Quantized:
        decode_ebchpf_quantized(request, result);
        break;
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
