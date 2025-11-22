#include "newcode/pipeline_runner.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <mutex>
#include <streambuf>

#include "newcode/awgn.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/decoder_api.hpp"
#include "newcode/info_extract.hpp"
#include "newcode/llr_known_prefix.hpp"
#include "newcode/llr_qpack.hpp"
#include "newcode/ofec_decoder.hpp"
#include "newcode/ofec_encoder.hpp"
#include "newcode/ofec_llr_matrix.hpp"
#include "newcode/qam.hpp"
#include "newcode/qam_llr.hpp"
#include "newcode/qfloat.hpp"
#include "newcode/interleaver.hpp"

namespace newcode {
namespace {

class NullStreamBuf : public std::streambuf {
 public:
  int overflow(int ch) override { return traits_type::not_eof(ch); }
};

std::ostream& null_stream() {
  static NullStreamBuf buf;
  static std::ostream stream(&buf);
  return stream;
}

void ensure_decoders_registered() {
  static std::once_flag once;
  std::call_once(once, [] {
    register_decoder_plain_factory();
    register_decoder_ebchPF_factory();
  });
}

std::vector<double> compute_early_stop_percentages(const std::vector<TileEarlyStopCounter>& counters)
{
  std::vector<double> pct;
  pct.reserve(counters.size());
  for (const auto& counter : counters) {
    double value = 0.0;
    if (counter.total > 0) {
      value = static_cast<double>(counter.triggered) /
              static_cast<double>(counter.total) * 100.0;
    }
    pct.push_back(value);
  }
  return pct;
}

// Flatten an encoded matrix (row-major) to a 0/1 bitstream.
std::vector<uint8_t> flatten_row_major(const Matrix<uint8_t>& matrix)
{
  std::vector<uint8_t> out;
  out.reserve(matrix.rows() * matrix.cols());
  for (size_t r = 0; r < matrix.rows(); ++r)
    for (size_t c = 0; c < matrix.cols(); ++c)
      out.push_back(matrix[r][c] & 1u);
  return out;
}

// 将硬比特矩阵(0/1)转换为“理想”LLR矩阵：0 -> +A，1 -> -A（供提取 TX 参考信息）
static Matrix<float> hard_bits_to_llr_matrix(const Matrix<uint8_t>& bits_mat, float A = 50.0f)
{
  Matrix<float> m(bits_mat.rows(), bits_mat.cols());
  for (size_t r = 0; r < bits_mat.rows(); ++r)
    for (size_t c = 0; c < bits_mat.cols(); ++c)
      m[r][c] = (bits_mat[r][c] ? -A : +A);
  return m;
}

struct LlrMode {
  LlrFormat format;
  std::size_t quant_bits;
};

LlrMode pick_llr_mode(const Params& params) {
  if (params.LLR_BITS == 16) {
    return {LlrFormat::Float, 16};
  }
  if (params.LLR_BITS >= 2 && params.LLR_BITS <= 15) {
    return {LlrFormat::Quantized, params.LLR_BITS};
  }
  throw std::runtime_error("[ERROR] Params::LLR_BITS must be 2..16");
}

} // namespace

PipelineResult run_pipeline(const Params& params,
                            const PipelineConfig& config,
                            const std::string& label,
                            float ebn0_dB)
{
  std::ostream& log = config.quiet ? null_stream() : std::cout;
  const bool verbose = !config.quiet;
  if (verbose) {
    log << "\n[RUN] Scenario: " << label << "\n";
  }
  
  // 生成信息比特
  auto info_bits = generate_bits(params);
  if (verbose) {
    log << "[INFO] (" << label << ") Generated bits: " << info_bits.size() << "\n";
  }

  // oFEC 编码
  auto code_matrix = ofec_encode(info_bits, params);
  if (verbose) {
    log << "[INFO] (" << label << ") oFEC matrix: " << code_matrix.rows()
        << " x " << code_matrix.cols() << "\n";
  }

  //抽取Tx发射参考信息
  const float TX_REF_LLR = 50.0f; // 任意足够大的幅度即可
  Matrix<float> tx_llr_mat = hard_bits_to_llr_matrix(code_matrix, TX_REF_LLR);
  auto tx_info_bits_ref = rx_info_from_bit_llr(tx_llr_mat, params);
  if (verbose) {
    log << "[INFO] (" << label << ") tx_info_bits_ref (by extractor): "
        << tx_info_bits_ref.size() << "\n";
  }

  // oFEC矩阵比特展平
  auto coded_bits = flatten_row_major(code_matrix);
  if (verbose) {
    log << "[INFO] (" << label << ") Coded bits (flattened): " << coded_bits.size() << "\n";
  }

  // 交织
  const std::size_t block_dim = static_cast<std::size_t>(Params::BITS_PER_SUBBLOCK_DIM);
  auto interleaver = newcode::Interleaver::build_from_shape(
      code_matrix.rows(), code_matrix.cols(), block_dim, block_dim, config.interleaver_name);

  auto coded_bits_itlv = interleaver.interleave_chunks(coded_bits);
  if (verbose) {
    log << "[INFO] (" << label << ") Interleaved bits: " << coded_bits_itlv.size() << "\n";
  }

  //调制
  unsigned n_bps = config.bits_per_symbol;
  if (n_bps == 0) n_bps = 2;
  std::string modulation_name;
  if (n_bps == 1) {
    modulation_name = "BPSK";
  } else if ((n_bps & 1u) == 0u) {
    modulation_name = (n_bps == 2)
        ? "QPSK"
        : (std::to_string(1ull << n_bps) + "-QAM");
  } else {
    throw std::invalid_argument("[ERROR] bits_per_symbol must be 1 or a positive even number.");
  }

  auto tx_syms = qam_modulate(coded_bits_itlv, n_bps);
  if (verbose) {
    log << "[INFO] (" << label << ") Modulation: " << modulation_name
        << " (n_bps=" << n_bps << ")\n";
    log << "[INFO] (" << label << ") Modulated symbols: " << tx_syms.size() << " (Es≈1)\n";
  }

  // 信道：加性高斯白噪声（AWGN）
  const int   N        = static_cast<int>(params.NUM_SUBBLOCK_COLS * params.BITS_PER_SUBBLOCK_DIM);
  const int   K        = 239;
  const int   TAKEBITS = K - N;
  const float code_rate = static_cast<float>(TAKEBITS) / static_cast<float>(N);
  const uint32_t awgn_seed = static_cast<uint32_t>(params.CHANNEL_SEED);
  auto rx_syms = add_awgn(tx_syms, ebn0_dB, n_bps, awgn_seed);

  if (verbose) {
    log << "[INFO] (" << label << ") Eb/N0 set to " << ebn0_dB << " dB\n";

    log << "[INFO] (" << label << ") Example symbols (TX -> RX):\n";
    for (size_t i = 0; i < std::min<size_t>(3, tx_syms.size()); ++i) {
      log << "  " << i
          << ": (" << tx_syms[i].real() << ", " << tx_syms[i].imag() << ")"
          << " -> (" << rx_syms[i].real() << ", " << rx_syms[i].imag() << ")\n";
    }
  }

  // qam解调计算信道 LLR
  auto llr = qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate);
  if (verbose) {
    log << "[INFO] (" << label << ") LLR count: " << llr.size()
        << " (should be tx_syms.size()*n_bps)\n";
    log << "[INFO] (" << label << ") First few LLRs: ";
    for (size_t i = 0; i < std::min<size_t>(8, llr.size()); ++i)
      log << llr[i] << (i + 1 < std::min<size_t>(8, llr.size()) ? ", " : "\n");
  }

  // 反交织
  auto llr_deint = interleaver.deinterleave_chunks(llr);

  // 转换为矩阵形式
  Matrix<float> llr_mat = llr_to_matrix_row_major(llr_deint, code_matrix.rows(), code_matrix.cols());

  // 已知前缀置零处理以及信道llr归一化
  apply_known_zero_prefix(llr_mat, params);

  // 构造解码器
  ensure_decoders_registered();
  auto decoder = make_decoder(config.decoder_name);
  if (!decoder) {
    throw std::runtime_error("[ERROR] make_decoder: unknown decoder '" + config.decoder_name + "'");
  }

  //确定量化所用的clip
  const auto llr_mode = pick_llr_mode(params);
  float quant_clip = 0.0f;
  if (llr_mode.format == LlrFormat::Quantized) {
    if (params.LLR_CLIP_RATIO > 0.0f) {
      std::vector<float> llr_values;
      llr_values.reserve(llr_mat.rows() * llr_mat.cols());
      for (size_t r = 0; r < llr_mat.rows(); ++r)
        for (size_t c = 0; c < llr_mat.cols(); ++c)
          llr_values.push_back(llr_mat[r][c]);
      quant_clip = compute_clip_from_ratio(llr_values.begin(), llr_values.end(),
                                           params.LLR_CLIP_RATIO);
    }
    if (quant_clip <= 0.0f) {
      quant_clip = params.LLR_CLIP;
    }
  }

  // 构造解码请求
  DecodeRequest request{
      .label = label,
      .channel_llr = llr_mat,
      .tx_llr_ref = &tx_llr_mat,
      .params = params,
      .format = llr_mode.format,
      .quant_bits = llr_mode.quant_bits,
      .quant_clip = quant_clip,
      .normalize_extrinsic = config.normalize_extrinsic,
      .quiet = config.quiet,
      .dump_quantized_llr = config.dump_quantized_llr,
      .quantized_llr_output_path = config.quantized_llr_output_path,
      .dump_float_llr = config.dump_quantized_llr,
      .float_llr_output_path = {},
      .dump_quantized_codes = config.dump_quantized_llr,
      .quantized_codes_output_path = {},
      .dump_work_llr = config.dump_work_llr,
      .work_llr_output_path = config.work_llr_output_path
  };

  // 执行解码
  auto decode_result = decoder->decode(request);

  // 提取解码后信息比特
  auto rx_info_bits_pre  = rx_info_from_bit_llr(decode_result.pre_decoder_llr,  params);
  auto rx_info_bits_post = rx_info_from_bit_llr(decode_result.post_decoder_llr, params);

  if (verbose) {
    log << "[INFO] (" << label << ") rx_info_bits: " << rx_info_bits_pre.size()
        << " (flattened, warmup skipped)\n";
  }

  const std::string pre_label  = label + " Pre-FEC";
  const std::string post_label = label + " Post-FEC";

  PipelineResult result;
  result.ebn0_db  = ebn0_dB;
  result.pre_fec_error_positions.clear();
  result.post_fec_error_positions.clear();
  result.pre_fec  = compute_and_print_ber(tx_info_bits_ref, rx_info_bits_pre,  pre_label.c_str(),
                                          params, &result.pre_fec_error_positions, config.quiet);
  result.post_fec = compute_and_print_ber(tx_info_bits_ref, rx_info_bits_post, post_label.c_str(),
                                          params, &result.post_fec_error_positions, config.quiet);
  result.tile_early_stop_pct = compute_early_stop_percentages(decode_result.tile_stats);
  result.dequantized_llr_path = decode_result.dequantized_llr_path;
  result.float_llr_path = decode_result.float_llr_path;
  result.quantized_codes_path = decode_result.quantized_codes_path;

  if (verbose) {
    log << "[DONE] (" << label << ") Pipeline bits -> channel -> decoder(" << config.decoder_name
        << ") using interleaver '" << config.interleaver_name << "' completed\n";
  }
  return result;
}

PipelineResult run_pipeline(const Params& params,
                            const std::string& label,
                            float ebn0_dB)
{
  PipelineConfig cfg{};
  return run_pipeline(params, cfg, label, ebn0_dB);
}

} // namespace newcode
