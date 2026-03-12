#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"
#include "newcode/io/dualwriter.hpp"
namespace ofec_single {

struct Config {
  std::string label;
  float ebn0_db;
  int chaseL_override;
  bool normalize_extrinsic;
  unsigned bits_per_symbol;
  int bitgen_seed;
  int channel_seed;
  bool enable_early_stop = true;
  int early_stop_detect_mode = 1;
  float alpha_fill;
  float beta_fill;
  std::vector<float> alpha_explicit;
  std::vector<float> beta_explicit;
  std::vector<int> siso_active_list;
  int mux_group_g = 1;
  bool mux_enable_reconfig = false;
  std::vector<newcode::mux::MuxEdge> mux_extra_bypass_edges;
  std::string interleaver_name;
  std::string decoder_name;
  bool generate_random_bits = true;
  bool normalize_known_prefix_tail = true;
  float quant_clip_ratio = 0.0f;
  std::size_t llr_bits = 16;
  bool dump_quantized_llr = false;
  std::string quantized_llr_output_path;
  bool dump_work_llr = false;
  std::string work_llr_output_path;
  newcode::Params::DebugTraceConfig debug_trace;
};



int run_ofec_single(const Config& config);

namespace detail {



std::optional<newcode::Params> build_params(const Config& cfg,
                                            io::DualWriter& log);

newcode::PipelineConfig build_pipeline_config(const Config& cfg);

void log_run_overview(const Config& cfg,
                      const newcode::Params& params,
                      io::DualWriter& log);

void log_pipeline_results(const newcode::PipelineResult& result,
                          io::DualWriter& log);

} 

}  
