#include "tpc_internal.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

#include "newcode/channel/awgn.hpp"
#include "newcode/ofec_llr_matrix.hpp"
#include "newcode/tx/bitgen/bitgen.hpp"
#include "newcode/tx/mod/qam.hpp"
#include "newcode/rx/demod/qam_llr.hpp"

namespace newcode {
namespace tpc {
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

uint8_t hard_decide(float value) {
  return (value >= 0.0f) ? 0u : 1u;
}

// TPC 信息位在左上 239x239，直接硬判决提取。
std::vector<uint8_t> extract_info_bits_from_llr(const matrix::Matrix<float>& llr) {
  if (llr.rows() < kTpcInfoDim || llr.cols() < kTpcInfoDim) {
    throw std::invalid_argument("extract_info_bits_from_llr: LLR matrix too small.");
  }
  std::vector<uint8_t> bits;
  bits.reserve(kTpcInfoDim * kTpcInfoDim);
  for (std::size_t r = 0; r < kTpcInfoDim; ++r) {
    for (std::size_t c = 0; c < kTpcInfoDim; ++c) {
      bits.push_back(hard_decide(llr[r][c]));
    }
  }
  return bits;
}

BerStats compute_ber_full(const std::vector<uint8_t>& ref_bits,
                          const std::vector<uint8_t>& rx_bits) {
  const std::size_t L = std::min(ref_bits.size(), rx_bits.size());
  BerStats stats;
  stats.total = L;
  for (std::size_t i = 0; i < L; ++i) {
    if ((ref_bits[i] ^ rx_bits[i]) & 1u) {
      ++stats.errors;
    }
  }
  stats.ber = (stats.total == 0) ? 0.0 : static_cast<double>(stats.errors) /
                                       static_cast<double>(stats.total);
  return stats;
}

// 强制生成 239x239 信息位，忽略 Params::NUM_INFO_BITS 的默认值。
std::vector<uint8_t> generate_info_bits(const Params& params) {
  Params local = params;
  local.NUM_INFO_BITS = kTpcInfoDim * kTpcInfoDim;
  return bitgen::generate_bits(local);
}

}  // namespace

TpcPipelineResult run_tpc_pipeline(const Params& params,
                                   const TpcPipelineConfig& config,
                                   const std::string& label,
                                   float ebn0_db) {
  if (params.LLR_BITS != 16) {
    throw std::invalid_argument("tpc_pipeline: only LLR_BITS=16 supported.");
  }
  if (config.num_blocks < 1) {
    throw std::invalid_argument("tpc_pipeline: num_blocks must be >= 1.");
  }
  std::ostream& log = config.quiet ? null_stream() : std::cout;
  const bool verbose = !config.quiet;

  if (verbose) {
    log << "\n[RUN] TPC Scenario: " << label << "\n";
    if (config.num_blocks > 1) {
      log << "[INFO] (" << label << ") blocks=" << config.num_blocks << "\n";
    }
  }

  std::size_t pre_errors_total = 0;
  std::size_t post_errors_total = 0;
  std::size_t total_bits = 0;
  int iters_total = 0;
  int early_stop_blocks = 0;
  int last_iters = 0;
  bool last_early_stop = false;
  std::vector<int> early_stop_hist(static_cast<std::size_t>(config.max_iters), 0);
  std::vector<double> per_block_ber_iter3;
  std::vector<double> per_block_ber_iter4;
  std::vector<int> per_block_early_stop_iter4;
  per_block_ber_iter3.reserve(static_cast<std::size_t>(config.num_blocks));
  per_block_ber_iter4.reserve(static_cast<std::size_t>(config.num_blocks));
  per_block_early_stop_iter4.reserve(static_cast<std::size_t>(config.num_blocks));
  const int iter3_index = 2;
  const int iter4_index = 3;
  const bool want_iter3 = config.max_iters > iter3_index;
  const bool want_iter4 = config.max_iters > iter4_index;

  for (int blk = 0; blk < config.num_blocks; ++blk) {
    Params local_params = params;
    local_params.BITGEN_SEED = params.BITGEN_SEED + blk;
    local_params.CHANNEL_SEED = params.CHANNEL_SEED + blk;
    const bool log_details = verbose && (blk == 0);

    // 生成信息比特并进行 TPC 行/列编码。
    const auto info_bits = generate_info_bits(local_params);
    if (log_details) {
      log << "[INFO] (" << label << ") Generated bits: " << info_bits.size() << "\n";
    }

    const auto code_matrix = tpc_encode(info_bits, local_params);
    if (log_details) {
      log << "[INFO] (" << label << ") TPC matrix: " << code_matrix.rows()
          << " x " << code_matrix.cols() << "\n";
    }

    // 码字矩阵按行展开后送调制与信道。
    auto coded_bits = matrix::flatten_row_major(
        code_matrix, [](uint8_t b) { return static_cast<uint8_t>(b & 1u); });
    if (log_details) {
      log << "[INFO] (" << label << ") Coded bits (flattened): " << coded_bits.size() << "\n";
    }

    unsigned n_bps = config.bits_per_symbol;
    if (n_bps == 0) {
      n_bps = 1;
    }
    if (n_bps != 1 && (n_bps & 1u) != 0u) {
      throw std::invalid_argument("tpc_pipeline: bits_per_symbol must be 1 or even.");
    }

    auto tx_syms = mod::qam_modulate(coded_bits, n_bps);
    if (log_details) {
      log << "[INFO] (" << label << ") Modulated symbols: " << tx_syms.size() << "\n";
    }

    const uint32_t awgn_seed = static_cast<uint32_t>(local_params.CHANNEL_SEED);
    auto rx_syms = channel::add_awgn(tx_syms, ebn0_db, n_bps, awgn_seed);

    // 码率为 (K/N)^2。
    const float rate = static_cast<float>(kTpcInfoDim) /
                       static_cast<float>(kTpcCodeDim);
    const float code_rate = rate * rate;

    auto llr = demod::qam_llr_from_ebn0(rx_syms, n_bps, ebn0_db, code_rate);
    if (log_details) {
      log << "[INFO] (" << label << ") LLR count: " << llr.size() << "\n";
    }

    matrix::Matrix<float> llr_mat =
        llr_to_matrix_row_major(llr, kTpcCodeDim, kTpcCodeDim);

    matrix::Matrix<float> iter3_llr;
    matrix::Matrix<float> iter4_llr;
    std::vector<bool> early_stop_start_flags(static_cast<std::size_t>(config.max_iters), false);
    TpcIterTrace trace;
    trace.report_iter_a = want_iter3 ? iter3_index : -1;
    trace.report_iter_b = want_iter4 ? iter4_index : -1;
    trace.iter_a_llr = want_iter3 ? &iter3_llr : nullptr;
    trace.iter_b_llr = want_iter4 ? &iter4_llr : nullptr;
    trace.early_stop_start_flags = &early_stop_start_flags;
    trace.continue_on_early_stop = true;

    const auto decode_result = tpc_decode_plain(
        llr_mat,
        local_params,
        config.max_iters,
        &config.alpha_schedule,
        &config.beta_schedule,
        &trace);

    const auto pre_bits = extract_info_bits_from_llr(decode_result.pre_decoder_llr);
    const auto post_bits = extract_info_bits_from_llr(decode_result.post_decoder_llr);

    const auto pre_stats = compute_ber_full(info_bits, pre_bits);
    const auto post_stats = compute_ber_full(info_bits, post_bits);
    double ber_iter3 = std::numeric_limits<double>::quiet_NaN();
    double ber_iter4 = std::numeric_limits<double>::quiet_NaN();
    if (want_iter3 && iter3_llr.rows() == kTpcCodeDim && iter3_llr.cols() == kTpcCodeDim) {
      const auto bits_iter3 = extract_info_bits_from_llr(iter3_llr);
      ber_iter3 = compute_ber_full(info_bits, bits_iter3).ber;
    }
    if (want_iter4 && iter4_llr.rows() == kTpcCodeDim && iter4_llr.cols() == kTpcCodeDim) {
      const auto bits_iter4 = extract_info_bits_from_llr(iter4_llr);
      ber_iter4 = compute_ber_full(info_bits, bits_iter4).ber;
    }
    per_block_ber_iter3.push_back(ber_iter3);
    per_block_ber_iter4.push_back(ber_iter4);
    int early_stop_iter4 = 0;
    if (want_iter4 &&
        early_stop_start_flags.size() > static_cast<std::size_t>(iter4_index)) {
      early_stop_iter4 = early_stop_start_flags[static_cast<std::size_t>(iter4_index)] ? 1 : 0;
    }
    per_block_early_stop_iter4.push_back(early_stop_iter4);

    pre_errors_total += pre_stats.errors;
    post_errors_total += post_stats.errors;
    total_bits += pre_stats.total;
    iters_total += decode_result.iters;
    if (decode_result.early_stop) {
      ++early_stop_blocks;
      if (decode_result.early_stop_iter >= 0 &&
          decode_result.early_stop_iter < config.max_iters) {
        ++early_stop_hist[static_cast<std::size_t>(decode_result.early_stop_iter)];
      }
    }
    last_iters = decode_result.iters;
    last_early_stop = decode_result.early_stop;
  }

  TpcPipelineResult result;
  result.ebn0_db = ebn0_db;
  result.pre_fec.errors = pre_errors_total;
  result.pre_fec.total = total_bits;
  result.pre_fec.ber = (total_bits == 0) ? 0.0
                                         : static_cast<double>(pre_errors_total) /
                                               static_cast<double>(total_bits);
  result.post_fec.errors = post_errors_total;
  result.post_fec.total = total_bits;
  result.post_fec.ber = (total_bits == 0) ? 0.0
                                          : static_cast<double>(post_errors_total) /
                                                static_cast<double>(total_bits);
  result.iters = last_iters;
  result.early_stop = (config.num_blocks == 1)
                          ? last_early_stop
                          : (early_stop_blocks == config.num_blocks);
  result.blocks = config.num_blocks;
  result.early_stop_blocks = early_stop_blocks;
  result.avg_iters = static_cast<double>(iters_total) /
                     static_cast<double>(config.num_blocks);
  result.early_stop_start_pct.clear();
  result.early_stop_start_pct.reserve(static_cast<std::size_t>(config.max_iters));
  for (int i = 0; i < config.max_iters; ++i) {
    const double pct = (config.num_blocks == 0)
                           ? 0.0
                           : static_cast<double>(early_stop_hist[static_cast<std::size_t>(i)]) *
                                 100.0 / static_cast<double>(config.num_blocks);
    result.early_stop_start_pct.push_back(pct);
  }
  result.per_block_ber_iter3 = std::move(per_block_ber_iter3);
  result.per_block_ber_iter4 = std::move(per_block_ber_iter4);
  result.per_block_early_stop_iter4 = std::move(per_block_early_stop_iter4);

  if (verbose) {
    log << "[RESULT] Pre-FEC BER=" << result.pre_fec.ber
        << " (errs=" << result.pre_fec.errors << "/" << result.pre_fec.total << ")\n";
    log << "[RESULT] Post-FEC BER=" << result.post_fec.ber
        << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")\n";
    log << "[INFO] (" << label << ") avg_iters=" << result.avg_iters
        << ", early_stop_blocks=" << result.early_stop_blocks
        << "/" << result.blocks << "\n";
    if (!result.early_stop_start_pct.empty()) {
      log << "[INFO] (" << label << ") early_stop_start_pct: ";
      for (std::size_t i = 0; i < result.early_stop_start_pct.size(); ++i) {
        log << result.early_stop_start_pct[i];
        if (i + 1 < result.early_stop_start_pct.size()) {
          log << ", ";
        }
      }
      log << "\n";
    }
  }

  return result;
}

}  // namespace tpc
}  // namespace newcode
