#pragma once

#include <cstddef>
#include <vector>

#include "newcode/common/matrix/matrix.hpp"
#include "newcode/params.hpp"
#include "newcode/tpc_pipeline.hpp"

namespace newcode {
namespace tpc {

inline constexpr std::size_t kTpcInfoDim = Params::BCH_K;
inline constexpr std::size_t kTpcCodeDim = Params::BCH_N;

matrix::Matrix<uint8_t> tpc_encode(const std::vector<uint8_t>& info_bits,
                                   const Params& params);

struct TpcDecodeResult {
  matrix::Matrix<float> pre_decoder_llr;
  matrix::Matrix<float> post_decoder_llr;
  int iters = 0;
  bool early_stop = false;
  int early_stop_iter = -1;
};

struct TpcIterTrace {
  int report_iter_a = -1;
  int report_iter_b = -1;
  matrix::Matrix<float>* iter_a_llr = nullptr;
  matrix::Matrix<float>* iter_b_llr = nullptr;
  std::vector<bool>* early_stop_start_flags = nullptr;
  bool continue_on_early_stop = false;
};

TpcDecodeResult tpc_decode_plain(const matrix::Matrix<float>& channel_llr,
                                 const Params& params,
                                 int max_iters,
                                 const std::vector<float>* alpha_schedule,
                                 const std::vector<float>* beta_schedule,
                                 TpcIterTrace* trace);

}  // namespace tpc
}  // namespace newcode
