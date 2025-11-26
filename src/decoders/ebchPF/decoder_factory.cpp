#include "newcode/decoder_api.hpp"

#include <iostream>
#include <mutex>
#include <stdexcept>

#include "newcode/llr_qpack.hpp"
#include "newcode/qfloat.hpp"
#include "newcode/quantized_llr_dump.hpp"

namespace newcode {
namespace {

template <int NBITS>
void decode_ebchPF_qfloat(const DecodeRequest& request, DecodeResult& result) {
  using Q = qfloat<NBITS>;

  const float clip = (request.quant_clip > 0.0f) ? request.quant_clip : request.params.LLR_CLIP;
  Q::set_clip(clip);

  if (!request.quiet) {
    std::cout << "[INFO] (" << request.label << ") ebchPF decoder running in qfloat<"
              << NBITS << ">\n";
  }
  auto quantized = quantize_matrix_to_qfloat<NBITS>(request.channel_llr, clip); //对channel_llr量化
  if (!request.quiet) {
    std::cout << "[INFO] (" << request.label << ") Quant clip=" << clip
              << " levels ±" << Q::Q() << "\n";
  }

  // pre_decoder_llr 一直使用原始通道LLR，以保证 pre-BER 基于未裁剪/未量化的数据
  result.pre_decoder_llr = request.channel_llr;
  auto decoded = ofec_decode_llr_ebchPF(quantized, request.params, &result.tile_stats,
                                        request.normalize_extrinsic, request.tx_llr_ref);//对量化后的LLR进行解码
  result.post_decoder_llr = dequantize_matrix_from_qfloat(decoded, clip); //对解码后的量化LLR矩阵反量化，得到post_decoder_llr

  if (request.dump_quantized_codes) {
    auto codes_mat = cast_matrix_from_qfloat(quantized);
    result.quantized_codes_path = dump_quantized_llr(codes_mat, request, request.dump_quantized_codes,request.quantized_codes_output_path, "_quantized_codes");
    result.dequantized_llr_path = dump_quantized_llr(result.pre_decoder_llr, request, request.dump_quantized_llr,request.quantized_llr_output_path, "_dequantized");
    result.float_llr_path = dump_quantized_llr(request.channel_llr, request, request.dump_float_llr,request.float_llr_output_path, "_float");
  }
}

void decode_ebchPF_quantized(const DecodeRequest& request, DecodeResult& result) {
  switch (request.quant_bits) {
    case 2:  decode_ebchPF_qfloat<2 >(request, result); break;
    case 3:  decode_ebchPF_qfloat<3 >(request, result); break;
    case 4:  decode_ebchPF_qfloat<4 >(request, result); break;
    case 5:  decode_ebchPF_qfloat<5 >(request, result); break;
    case 6:  decode_ebchPF_qfloat<6 >(request, result); break;
    case 7:  decode_ebchPF_qfloat<7 >(request, result); break;
    case 8:  decode_ebchPF_qfloat<8 >(request, result); break;
    case 9:  decode_ebchPF_qfloat<9 >(request, result); break;
    case 10: decode_ebchPF_qfloat<10>(request, result); break;
    case 11: decode_ebchPF_qfloat<11>(request, result); break;
    case 12: decode_ebchPF_qfloat<12>(request, result); break;
    case 13: decode_ebchPF_qfloat<13>(request, result); break;
    case 14: decode_ebchPF_qfloat<14>(request, result); break;
    case 15: decode_ebchPF_qfloat<15>(request, result); break;
    default:
      throw std::runtime_error("[ERROR] Unsupported quant_bits for ebchPF decoder: " +
                               std::to_string(request.quant_bits));
  }
}

class ebchPFDecoder final : public IDecoder {
public:
  DecodeResult decode(const DecodeRequest& request) override {
    DecodeResult result;
    result.stats = {};

    switch (request.format) {
      case LlrFormat::Float:
        if (!request.quiet) {
          std::cout << "[INFO] (" << request.label << ") ebchPF decoder running in FLOAT\n";
        }
        if (!request.quiet) {
          if (request.dump_quantized_llr) {
            std::cout << "[WARN] (" << request.label << ") dump_quantized_llr ignored in FLOAT mode\n";
          }
          if (request.dump_float_llr) {
            std::cout << "[WARN] (" << request.label << ") dump_float_llr ignored in FLOAT mode (already float)\n";
          }
          if (request.dump_quantized_codes) {
            std::cout << "[WARN] (" << request.label << ") dump_quantized_codes ignored in FLOAT mode\n";
          }
        }
        result.pre_decoder_llr = request.channel_llr;
        result.post_decoder_llr =
            ofec_decode_llr_ebchPF(result.pre_decoder_llr, request.params, &result.tile_stats,
                                   request.normalize_extrinsic, request.tx_llr_ref);
        break;

      case LlrFormat::Quantized:
        decode_ebchPF_quantized(request, result);
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
        return std::make_unique<ebchPFDecoder>();
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
