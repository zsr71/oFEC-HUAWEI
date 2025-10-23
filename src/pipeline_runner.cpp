#include "newcode/pipeline_runner.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "newcode/awgn.hpp"
#include "newcode/bitgen.hpp"
#include "newcode/info_extract.hpp"
#include "newcode/llr_known_prefix.hpp"
#include "newcode/llr_qpack.hpp"
#include "newcode/ofec_decoder.hpp"
#include "newcode/ofec_encoder.hpp"
#include "newcode/ofec_llr_matrix.hpp"
#include "newcode/qam.hpp"
#include "newcode/qam_llr.hpp"
#include "newcode/qfloat.hpp"

namespace newcode {
namespace {

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

// Run the pipeline using plain float LLRs.
PipelineResult run_with_float(const Matrix<float>& llr_mat,
                              const Params& params,
                              const std::vector<uint8_t>& info_bits,
                              const std::string& label)
{
  std::cout << "[INFO] (" << label << ") Running decoder in FLOAT (no quantization, no clipping)\n";

  std::vector<TileEarlyStopCounter> tile_stats;
  Matrix<float> post_f = ofec_decode_llr(llr_mat, params, &tile_stats);

  auto rx_info_bits_pre  = rx_info_from_bit_llr(llr_mat, params);
  auto rx_info_bits_post = rx_info_from_bit_llr(post_f,  params);

  std::cout << "[INFO] (" << label << ") rx_info_bits: " << rx_info_bits_pre.size()
            << " (flattened, warmup skipped)\n";

  const std::string pre_label  = label + " Pre-FEC";
  const std::string post_label = label + " Post-FEC";

  PipelineResult result;
  result.pre_fec  = compute_and_print_ber(info_bits, rx_info_bits_pre,  pre_label.c_str(),  params);
  result.post_fec = compute_and_print_ber(info_bits, rx_info_bits_post, post_label.c_str(), params);
  result.tile_early_stop_pct = compute_early_stop_percentages(tile_stats);

  std::cout << "[DONE] (" << label << ") Pipeline(float) bits -> oFEC -> QAM -> AWGN -> QAM LLR -> decode(float) -> info extract & BER\n";
  return result;
}

// Run the pipeline using qfloat-quantised LLRs.
template<int NBITS>
PipelineResult run_with_qfloat(const Matrix<float>& llr_mat,
                               const Params& params,
                               const std::vector<uint8_t>& info_bits,
                               const std::string& label)
{
  auto q_mat = quantize_matrix_to_qfloat<NBITS>(llr_mat, qfloat<NBITS>::DEFAULT_CLIP);
  std::cout << "[INFO] (" << label << ") Quantized qfloat<" << NBITS << ">: clip="
            << qfloat<NBITS>::DEFAULT_CLIP
            << "  code range [-" << qfloat<NBITS>::Q() << ", +" << qfloat<NBITS>::Q() << "]\n";

  std::vector<TileEarlyStopCounter> tile_stats;
  auto q_dec = ofec_decode_llr(q_mat, params, &tile_stats);

  auto pre_f  = cast_matrix_from_qfloat(q_mat);
  auto post_f = cast_matrix_from_qfloat(q_dec);

  auto rx_info_bits_pre  = rx_info_from_bit_llr(pre_f,  params);
  auto rx_info_bits_post = rx_info_from_bit_llr(post_f, params);

  std::cout << "[INFO] (" << label << ") rx_info_bits: " << rx_info_bits_pre.size()
            << " (flattened, warmup skipped)\n";

  const std::string pre_label  = label + " Pre-FEC";
  const std::string post_label = label + " Post-FEC";

  PipelineResult result;
  result.pre_fec  = compute_and_print_ber(info_bits, rx_info_bits_pre,  pre_label.c_str(),  params);
  result.post_fec = compute_and_print_ber(info_bits, rx_info_bits_post, post_label.c_str(), params);
  result.tile_early_stop_pct = compute_early_stop_percentages(tile_stats);

  std::cout << "[DONE] (" << label << ") Pipeline(qfloat<" << NBITS
            << ">) bits -> oFEC -> QAM -> AWGN -> QAM LLR -> quant(qfloat) -> decode -> info extract & BER\n";
  return result;
}

} // namespace

PipelineResult run_pipeline(const Params& params, const std::string& label)
{
  std::cout << "\n[RUN] Scenario: " << label << "\n";

  auto info_bits = generate_bits(params);
  std::cout << "[INFO] (" << label << ") Generated bits: " << info_bits.size() << "\n";

  auto code_matrix = ofec_encode(info_bits, params);
  std::cout << "[INFO] (" << label << ") oFEC matrix: " << code_matrix.rows()
            << " x " << code_matrix.cols() << "\n";

  auto coded_bits = flatten_row_major(code_matrix);
  std::cout << "[INFO] (" << label << ") Coded bits (flattened): " << coded_bits.size() << "\n";

  const unsigned n_bps = 2; // QPSK
  auto tx_syms = qam_modulate(coded_bits, n_bps);
  std::cout << "[INFO] (" << label << ") Modulated symbols: " << tx_syms.size() << " (Es≈1)\n";

  const float ebn0_dB  = 3.24f;
  const int   N        = static_cast<int>(params.NUM_SUBBLOCK_COLS * params.BITS_PER_SUBBLOCK_DIM);
  const int   K        = 239;
  const int   TAKEBITS = K - N;
  const float code_rate = static_cast<float>(TAKEBITS) / static_cast<float>(N);
  const uint32_t awgn_seed = static_cast<uint32_t>(params.BITGEN_SEED + 100);
  auto rx_syms = add_awgn(tx_syms, ebn0_dB, n_bps, awgn_seed);

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

  Matrix<float> llr_mat = llr_to_matrix_row_major(llr, code_matrix.rows(), code_matrix.cols());

  apply_known_zero_prefix(llr_mat, params);

  if (params.LLR_BITS == 16) {
    return run_with_float(llr_mat, params, info_bits, label);
  } else if (params.LLR_BITS == 5) {
    return run_with_qfloat<5>(llr_mat, params, info_bits, label);
  } else if (params.LLR_BITS == 4) {
    return run_with_qfloat<4>(llr_mat, params, info_bits, label);
  }

  throw std::runtime_error("[ERROR] Unsupported p.LLR_BITS value");
}

} // namespace newcode
