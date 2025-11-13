#pragma once
#pragma once

#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "newcode/ofec_llr_matrix.hpp"
#include "newcode/ofec_decoder.hpp"
#include "newcode/params.hpp"

namespace newcode {

enum class LlrFormat {
  Float,
  Quantized
};

struct DecodeStats {
  int iters = 0;
  bool success = false;
};

struct DecodeRequest {
  std::string_view label;
  const Matrix<float>& channel_llr;
  const Params& params;
  LlrFormat format = LlrFormat::Float;
  std::size_t quant_bits = 16;
  float quant_clip = 0.0f;
  bool normalize_extrinsic = true;
  bool quiet = false;
};

struct DecodeResult {
  Matrix<float> pre_decoder_llr;
  Matrix<float> post_decoder_llr;
  std::vector<TileEarlyStopCounter> tile_stats;
  DecodeStats stats;
};

struct IDecoder {
  virtual ~IDecoder() = default;
  virtual DecodeResult decode(const DecodeRequest& request) = 0;
};

using DecoderFactory = std::function<std::unique_ptr<IDecoder>(const std::string& name)>;

void register_decoder_factory(DecoderFactory factory);

std::unique_ptr<IDecoder> make_decoder(const std::string& name);

// Built-in decoder registration hooks (implemented by each decoder module)
void register_decoder_plain_factory();
void register_decoder_ebchPF_factory();

} // namespace newcode
