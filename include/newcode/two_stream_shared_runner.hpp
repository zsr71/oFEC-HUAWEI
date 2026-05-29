#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "newcode/params.hpp"
#include "newcode/pipeline_runner.hpp"

namespace newcode::two_stream_shared {

struct StreamConfig {
  std::string label;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = DEFAULT_EBN0_DB;
};

struct Config {
  Params params;
  PipelineConfig pipeline;
  StreamConfig stream_a;
  StreamConfig stream_b;
};

struct StreamQuantizationStats {
  float quant_clip = 0.0f;
  std::size_t saturated_values = 0;
  std::size_t total_values = 0;
  double saturation_ratio = 0.0;
};

struct SharedCoreAggregateStats {
  std::size_t total_batches = 0;
  std::size_t total_rows = 0;
  std::size_t produced_rows = 0;
  std::size_t failed_rows = 0;
  std::vector<std::size_t> produced_rows_per_stream;
  std::vector<std::size_t> failed_rows_per_stream;
};

struct Observability {
  float shared_quant_clip = 0.0f;
  StreamQuantizationStats stream_a_quantization;
  StreamQuantizationStats stream_b_quantization;
  SharedCoreAggregateStats shared_core;
};

struct Result {
  PipelineResult stream_a;
  PipelineResult stream_b;
  Observability observability;
};

Result run_two_stream_shared(const Config& config);

}  // namespace newcode::two_stream_shared
