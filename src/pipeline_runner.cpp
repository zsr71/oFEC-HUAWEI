#include "newcode/pipeline_runner.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#include <mutex>

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
#include "newcode/interleaver.hpp"

namespace newcode {
namespace {

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

LlrFormat pick_llr_format(const Params& params) {
  if (params.LLR_BITS == 16) return LlrFormat::Float;
  if (params.LLR_BITS == 5)  return LlrFormat::QFloat5;
  if (params.LLR_BITS == 4)  return LlrFormat::QFloat4;
  throw std::runtime_error("[ERROR] Unsupported Params::LLR_BITS value");
}

} // namespace

PipelineResult run_pipeline(const Params& params,
                            const PipelineConfig& config,
                            const std::string& label,
                            float ebn0_dB)
{
  std::cout << "\n[RUN] Scenario: " << label << "\n";

  auto info_bits = generate_bits(params);
  std::cout << "[INFO] (" << label << ") Generated bits: " << info_bits.size() << "\n";

  // 编码
  auto code_matrix = ofec_encode(info_bits, params);
  std::cout << "[INFO] (" << label << ") oFEC matrix: " << code_matrix.rows()
            << " x " << code_matrix.cols() << "\n";

  // === 新增：从编码矩阵得到“理想”LLR，再用同一提取器抽取 TX 参考信息 ===
  const float TX_REF_LLR = 50.0f; // 任意足够大的幅度即可
  Matrix<float> tx_llr_mat = hard_bits_to_llr_matrix(code_matrix, TX_REF_LLR);
  auto tx_info_bits_ref = rx_info_from_bit_llr(tx_llr_mat, params);
  std::cout << "[INFO] (" << label << ") tx_info_bits_ref (by extractor): "
            << tx_info_bits_ref.size() << "\n";

  // 展平比特 -> 调制
  auto coded_bits = flatten_row_major(code_matrix);
  std::cout << "[INFO] (" << label << ") Coded bits (flattened): " << coded_bits.size() << "\n";

  const std::size_t block_dim = static_cast<std::size_t>(Params::BITS_PER_SUBBLOCK_DIM);
  auto interleaver = newcode::Interleaver::build_from_shape(
      code_matrix.rows(), code_matrix.cols(), block_dim, block_dim, config.interleaver_name);

  auto coded_bits_itlv = interleaver.interleave_chunks(coded_bits);
  std::cout << "[INFO] (" << label << ") Interleaved bits: " << coded_bits_itlv.size() << "\n";

  const unsigned n_bps = 2; // QPSK
  auto tx_syms = qam_modulate(coded_bits_itlv, n_bps);
  std::cout << "[INFO] (" << label << ") Modulated symbols: " << tx_syms.size() << " (Es≈1)\n";

  const int   N        = static_cast<int>(params.NUM_SUBBLOCK_COLS * params.BITS_PER_SUBBLOCK_DIM);
  const int   K        = 239;
  const int   TAKEBITS = K - N;
  const float code_rate = static_cast<float>(TAKEBITS) / static_cast<float>(N);
  const uint32_t awgn_seed = static_cast<uint32_t>(params.BITGEN_SEED + 100);
  auto rx_syms = add_awgn(tx_syms, ebn0_dB, n_bps, awgn_seed);

  std::cout << "[INFO] (" << label << ") Eb/N0 set to " << ebn0_dB << " dB\n";

  std::cout << "[INFO] (" << label << ") Example symbols (TX -> RX):\n";
  for (size_t i = 0; i < std::min<size_t>(3, tx_syms.size()); ++i) {
    std::cout << "  " << i
              << ": (" << tx_syms[i].real() << ", " << tx_syms[i].imag() << ")"
              << " -> (" << rx_syms[i].real() << ", " << rx_syms[i].imag() << ")\n";
  }

  auto llr = qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate);
  std::cout << "[INFO] (" << label << ") LLR count: " << llr.size()
            << " (should be tx_syms.size()*n_bps)\n";
  std::cout << "[INFO] (" << label << ") First few LLRs: ";
  for (size_t i = 0; i < std::min<size_t>(8, llr.size()); ++i)
    std::cout << llr[i] << (i + 1 < std::min<size_t>(8, llr.size()) ? ", " : "\n");

  auto llr_deint = interleaver.deinterleave_chunks(llr);

  Matrix<float> llr_mat = llr_to_matrix_row_major(llr_deint, code_matrix.rows(), code_matrix.cols());

  // 已知前缀的先验处理（如果使用）
  apply_known_zero_prefix(llr_mat, params);

  ensure_decoders_registered();
  auto decoder = make_decoder(config.decoder_name);
  if (!decoder) {
    throw std::runtime_error("[ERROR] make_decoder: unknown decoder '" + config.decoder_name + "'");
  }

  DecodeRequest request{
      .label = label,
      .channel_llr = llr_mat,
      .params = params,
      .format = pick_llr_format(params),
      .normalize_extrinsic = config.normalize_extrinsic
  };

  auto decode_result = decoder->decode(request);

  auto rx_info_bits_pre  = rx_info_from_bit_llr(decode_result.pre_decoder_llr,  params);
  auto rx_info_bits_post = rx_info_from_bit_llr(decode_result.post_decoder_llr, params);

  std::cout << "[INFO] (" << label << ") rx_info_bits: " << rx_info_bits_pre.size()
            << " (flattened, warmup skipped)\n";

  const std::string pre_label  = label + " Pre-FEC";
  const std::string post_label = label + " Post-FEC";

  PipelineResult result;
  result.ebn0_db  = ebn0_dB;
  result.pre_fec_error_positions.clear();
  result.post_fec_error_positions.clear();
  result.pre_fec  = compute_and_print_ber(tx_info_bits_ref, rx_info_bits_pre,  pre_label.c_str(),  params, &result.pre_fec_error_positions);
  result.post_fec = compute_and_print_ber(tx_info_bits_ref, rx_info_bits_post, post_label.c_str(), params, &result.post_fec_error_positions);
  result.tile_early_stop_pct = compute_early_stop_percentages(decode_result.tile_stats);

  std::cout << "[DONE] (" << label << ") Pipeline bits -> channel -> decoder(" << config.decoder_name
            << ") using interleaver '" << config.interleaver_name << "' completed\n";
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
