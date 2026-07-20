#include "tpc_internal.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>

#include "newcode/common/bch/bch_255_239.hpp"
#include "newcode/rx/ofec/chase/chase256.hpp"

namespace newcode {
namespace tpc {
namespace {

uint8_t hard_decide(float value) {
  return (value >= 0.0f) ? 0u : 1u;
}

bool check_bch_256(const std::array<uint8_t, kTpcCodeDim>& bits) {
  if (!bch::bch_255_239_syndromes_zero_cw_255(bits.data())) {
    return false;
  }
  uint8_t parity = 0;
  for (std::size_t i = 0; i < kTpcCodeDim - 1; ++i) {
    parity ^= bits[i];
  }
  return parity == bits[kTpcCodeDim - 1];
}

// 行列全部通过 BCH 综合征 + overall parity 才算早停。
bool all_syndrome_ok(const matrix::Matrix<float>& llr) {
  const std::size_t rows = llr.rows();
  const std::size_t cols = llr.cols();
  if (rows != kTpcCodeDim || cols != kTpcCodeDim) {
    throw std::invalid_argument("tpc_decode_plain: LLR matrix size mismatch.");
  }

  std::array<uint8_t, kTpcCodeDim> bits{};
  for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
    for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
      bits[c] = hard_decide(llr[r][c]);
    }
    if (!check_bch_256(bits)) {
      return false;
    }
  }
  for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
    for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
      bits[r] = hard_decide(llr[r][c]);
    }
    if (!check_bch_256(bits)) {
      return false;
    }
  }
  return true;
}

float pick_schedule_value(const std::vector<float>& schedule, int index) {
  return schedule[static_cast<std::size_t>(index)];
}

// 行方向 SISO Chase：输出外信息，按 y_next = Lch + alpha * omega 更新。
void decode_rows_plain(const matrix::Matrix<float>& input,
                       const matrix::Matrix<float>& channel,
                       matrix::Matrix<float>& output,
                       const Params& params) {
  const float alpha = params.ALPHA;
  std::array<float, kTpcCodeDim> lin{};
  std::array<float, kTpcCodeDim> omega{};
  for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
    for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
      lin[c] = input[r][c];
    }
    chase::chase_decode_256_plain<float>(lin.data(), omega.data(), params);
    for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
      output[r][c] = channel[r][c] + alpha * omega[c];
    }
  }
}

// 列方向 SISO Chase：输出外信息，按 y_next = Lch + alpha * omega 更新。
void decode_cols_plain(const matrix::Matrix<float>& input,
                       const matrix::Matrix<float>& channel,
                       matrix::Matrix<float>& output,
                       const Params& params) {
  const float alpha = params.ALPHA;
  std::array<float, kTpcCodeDim> lin{};
  std::array<float, kTpcCodeDim> omega{};
  for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
    for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
      lin[r] = input[r][c];
    }
    chase::chase_decode_256_plain<float>(lin.data(), omega.data(), params);
    for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
      output[r][c] = channel[r][c] + alpha * omega[r];
    }
  }
}

void decode_cols_plain_with_base(const matrix::Matrix<float>& input,
                                 const matrix::Matrix<float>& base,
                                 matrix::Matrix<float>& output,
                                 const Params& params) {
  const float alpha = params.ALPHA;
  std::array<float, kTpcCodeDim> lin{};
  std::array<float, kTpcCodeDim> omega{};
  for (std::size_t c = 0; c < kTpcCodeDim; ++c) {
    for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
      lin[r] = input[r][c];
    }
    chase::chase_decode_256_plain<float>(lin.data(), omega.data(), params);
    for (std::size_t r = 0; r < kTpcCodeDim; ++r) {
      output[r][c] = base[r][c] + alpha * omega[r];
    }
  }
}

}  // namespace

TpcDecodeResult tpc_decode_plain(const matrix::Matrix<float>& channel_llr,
                                 const Params& params,
                                 int max_iters,
                                 const std::vector<float>* alpha_schedule,
                                 const std::vector<float>* beta_schedule,
                                 TpcIterTrace* trace) {
  if (channel_llr.rows() != kTpcCodeDim || channel_llr.cols() != kTpcCodeDim) {
    throw std::invalid_argument("tpc_decode_plain: LLR matrix size mismatch.");
  }
  if (max_iters < 1) {
    throw std::invalid_argument("tpc_decode_plain: max_iters must be >= 1.");
  }
  if (!alpha_schedule || !beta_schedule) {
    throw std::invalid_argument("tpc_decode_plain: schedules must be provided.");
  }
  const int schedule_len = max_iters * 2;
  if (static_cast<int>(alpha_schedule->size()) != schedule_len) {
    throw std::invalid_argument("tpc_decode_plain: alpha_schedule size mismatch.");
  }
  if (static_cast<int>(beta_schedule->size()) != schedule_len) {
    throw std::invalid_argument("tpc_decode_plain: beta_schedule size mismatch.");
  }

  TpcDecodeResult result;
  result.pre_decoder_llr = channel_llr;

  matrix::Matrix<float> work = channel_llr;
  matrix::Matrix<float> scratch(kTpcCodeDim, kTpcCodeDim);
  Params step_params = params;

  // 行/列交替迭代，综合征全通过即提前停止。
  for (int iter = 0; iter < max_iters; ++iter) {
    const bool start_ok = all_syndrome_ok(work);
    if (trace && trace->early_stop_start_flags &&
        trace->early_stop_start_flags->size() > static_cast<std::size_t>(iter)) {
      (*trace->early_stop_start_flags)[static_cast<std::size_t>(iter)] = start_ok;
    }
    if (start_ok && !result.early_stop) {
      result.early_stop = true;
      result.early_stop_iter = iter;
    }
    if (start_ok && (!trace || !trace->continue_on_early_stop)) {
      result.iters = iter;
      //break;
    }
    const int row_step = iter * 2;
    const int col_step = row_step + 1;
    step_params.ALPHA = pick_schedule_value(*alpha_schedule, row_step);
    step_params.beta = pick_schedule_value(*beta_schedule, row_step);
    decode_rows_plain(work, channel_llr, scratch, step_params);
    step_params.ALPHA = pick_schedule_value(*alpha_schedule, col_step);
    step_params.beta = pick_schedule_value(*beta_schedule, col_step);
    if (iter == max_iters - 1) {
      decode_cols_plain_with_base(scratch, scratch, work, step_params);
    } else {
      decode_cols_plain(scratch, channel_llr, work, step_params);
    }
    result.iters = iter + 1;
    if (trace) {
      if (trace->iter_a_llr && iter == trace->report_iter_a) {
        *trace->iter_a_llr = work;
      }
      if (trace->iter_b_llr && iter == trace->report_iter_b) {
        *trace->iter_b_llr = work;
      }
    }
  }

  result.post_decoder_llr = work;
  return result;
}

}  // namespace tpc
}  // namespace newcode
