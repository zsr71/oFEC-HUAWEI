#pragma once

#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"

namespace ofec_single {

struct Config {
  std::string label;
  float ebn0_db;
  int chaseL_override;
  bool normalize_extrinsic;
  unsigned bits_per_symbol;
  int bitgen_seed;
  int channel_seed;
  float alpha_fill;
  float beta_fill;
  std::vector<float> alpha_explicit;
  std::vector<float> beta_explicit;
  std::string interleaver_name;
  std::string decoder_name;
  bool generate_random_bits = true;
  bool normalize_known_prefix_tail = true;
  float quant_clip_ratio = 0.0f;
  std::size_t llr_bits = 16;
  newcode::Params::DebugTraceConfig debug_trace;
};

class DualWriter {
 public:
  explicit DualWriter(std::ofstream& file);

  template <typename T>
  DualWriter& operator<<(const T& value) {
    if (console_) {
      *console_ << value;
    }
    if (file_ && file_->is_open()) {
      *file_ << value;
    }
    return *this;
  }

  DualWriter& operator<<(std::ostream& (*manip)(std::ostream&));
  DualWriter& operator<<(std::ios_base& (*manip)(std::ios_base&));

 private:
  std::ostream* console_;
  std::ofstream* file_;
};

int run_ofec_single(const Config& config);

namespace detail {

std::ofstream prepare_log_file(const std::filesystem::path& data_dir,
                               std::string& log_path);

std::optional<newcode::Params> build_params(const Config& cfg,
                                            DualWriter& log);

newcode::PipelineConfig build_pipeline_config(const Config& cfg);

void log_run_overview(const Config& cfg,
                      const newcode::Params& params,
                      DualWriter& log);

void log_pipeline_results(const newcode::PipelineResult& result,
                          DualWriter& log);

}  // namespace detail

}  // namespace ofec_single
