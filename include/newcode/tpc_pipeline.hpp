#pragma once

#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/rx/ber/ber.hpp"

namespace newcode {
namespace tpc {

struct TpcPipelineConfig {
  unsigned bits_per_symbol = 1;
  int max_iters = 4;
  int num_blocks = 1;
  std::vector<float> alpha_schedule;
  std::vector<float> beta_schedule;
  bool quiet = false;
};

struct TpcPipelineResult {
  float ebn0_db = 0.0f;
  BerStats pre_fec;
  BerStats post_fec;
  int iters = 0;
  bool early_stop = false;
  int blocks = 1;
  int early_stop_blocks = 0;
  double avg_iters = 0.0;
  std::vector<double> early_stop_start_pct;
  std::vector<double> per_block_ber_iter3;
  std::vector<double> per_block_ber_iter4;
  std::vector<int> per_block_early_stop_iter4;
};

TpcPipelineResult run_tpc_pipeline(const Params& params,
                                   const TpcPipelineConfig& config,
                                   const std::string& label,
                                   float ebn0_db);

}  // namespace tpc
}  // namespace newcode
