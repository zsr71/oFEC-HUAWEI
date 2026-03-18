#include "new_float_only/pipeline_runner.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

#include "new_float_only/channel/awgn.hpp"
#include "new_float_only/common/matrix/hard_bits_to_llr_matrix.hpp"
#include "new_float_only/common/matrix/info_extract.hpp"
#include "new_float_only/common/matrix/matrix.hpp"
#include "new_float_only/llr_known_prefix.hpp"
#include "new_float_only/ofec_decoder.hpp"
#include "new_float_only/ofec_llr_matrix.hpp"
#include "new_float_only/rx/demod/qam_llr.hpp"
#include "new_float_only/tx/bitgen/bitgen.hpp"
#include "new_float_only/tx/mod/qam.hpp"
#include "new_float_only/tx/ofecencoder/ofec_encoder.hpp"

namespace new_float_only {
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

}  // namespace

/**
 * 执行一次完整的发端/收端浮点链路。
 * 关键步骤依次为：生成原始比特、oFEC 编码、调制、加噪、解调、矩阵化、
 * 前缀处理、plain 浮点解码、信息位提取、BER 统计。
 * 当前子工程固定为“无交织直通”链路，因此 coded bits 和解调 LLR 都不会再经过重排。
 */
PipelineResult run_pipeline(const Params& params,
                            const PipelineConfig& config,
                            const std::string& label,
                            float ebn0_dB) {
  std::ostream& log = config.quiet ? null_stream() : std::cout;
  const bool verbose = !config.quiet;

  if (verbose) {
    log << "\n[RUN] Scenario: " << label << "\n";
  }

  const std::vector<uint8_t> info_bits =
      bitgen::generate_bits(params.NUM_INFO_BITS,
                            config.bitgen_seed,
                            config.generate_random_bits);

  const auto code_matrix = ofecencoder::ofec_encode(info_bits, params);

  // 发送端参考序列仍使用当前提取器定义，便于和旧 plain+float 统计口径对齐。
  const float tx_ref_llr = 50.0f;
  const matrix::Matrix<float> tx_llr_mat =
      hard_bits_to_llr_matrix(code_matrix, tx_ref_llr);
  const auto tx_info_bits_ref = matrix::rx_info_from_bit_llr(tx_llr_mat, params);

  const auto coded_bits =
      flatten_row_major(code_matrix, [](uint8_t bit) { return static_cast<uint8_t>(bit & 1u); });

  unsigned n_bps = config.bits_per_symbol;
  if (n_bps == 0) {
    n_bps = 2;
  }
  if (n_bps != 1 && (n_bps & 1u)) {
    throw std::invalid_argument("bits_per_symbol must be 1 or a positive even number");
  }

  // 发端直接按编码矩阵的行优先顺序送入调制器，不做额外交织。
  const auto tx_syms = mod::qam_modulate(coded_bits, n_bps);

  const int n = static_cast<int>(Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM);
  const int take_bits = static_cast<int>(Params::BCH_K) - n;
  const float code_rate = static_cast<float>(take_bits) / static_cast<float>(n);
  const auto rx_syms =
      channel::add_awgn(tx_syms, ebn0_dB, n_bps, static_cast<uint32_t>(config.channel_seed));

  const auto llr = demod::qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate);

  // 收端也按解调输出的原顺序直接回填矩阵，不再做反交织。
  matrix::Matrix<float> llr_mat =
      llr_to_matrix_row_major(llr, code_matrix.rows(), code_matrix.cols());

  apply_known_zero_prefix(llr_mat, params);

  // 这里直接调用新的 float-only 解码入口，不再经过运行时 decoder registry。
  const matrix::Matrix<float> post_decoder_llr =
      decode_plain_llr(llr_mat, params, config.normalize_extrinsic, &tx_llr_mat);

  const auto rx_info_bits_pre = matrix::rx_info_from_bit_llr(llr_mat, params);
  const auto rx_info_bits_post = matrix::rx_info_from_bit_llr(post_decoder_llr, params);

  PipelineResult result;
  result.ebn0_db = ebn0_dB;
  std::vector<std::size_t>* pre_error_positions =
      config.collect_error_positions ? &result.pre_fec_error_positions : nullptr;
  std::vector<std::size_t>* post_error_positions =
      config.collect_error_positions ? &result.post_fec_error_positions : nullptr;
  result.pre_fec = compute_and_print_ber(tx_info_bits_ref,
                                         rx_info_bits_pre,
                                         (label + " Pre-FEC").c_str(),
                                         params,
                                         pre_error_positions,
                                         config.quiet);
  result.post_fec = compute_and_print_ber(tx_info_bits_ref,
                                          rx_info_bits_post,
                                          (label + " Post-FEC").c_str(),
                                          params,
                                          post_error_positions,
                                          config.quiet);
  if (params.DUMP_WORK_LLR) {
    result.work_llr_path = params.WORK_LLR_OUTPUT_PATH;
  }

  if (verbose) {
    log << "[DONE] (" << label << ") Float plain pipeline finished\n";
  }
  return result;
}

}  // namespace new_float_only
