#pragma once

#include <cstddef>
#include <cstdint>
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

enum class SharedHybridClass : std::uint8_t {
  None = 0,
  BchHardDecoded,
  Clean,
  ParityOnly,
  OneMain,
  OneMainPlusParity,
  TwoMain,
  Suspicious,
  HardFail
};

enum class SharedRowFinalTag : std::uint8_t {
  SoftDecode = 0,
  EarlyStopAction,
  HardFinish,
  Unscheduled
};

struct SharedTileSample {
  std::size_t invocation = 0;
  std::size_t tile_index = 0;
  std::size_t stream_rows_a = 0;
  std::size_t stream_rows_b = 0;
  std::size_t rows_total = 0;
  std::size_t rows_early_stop = 0;
  std::size_t rows_not_early_stop = 0;
  std::size_t rows_hard_finish = 0;
  std::size_t rows_need_siso_before_mux = 0;
  std::size_t rows_soft_scheduled = 0;
  std::size_t rows_unscheduled = 0;
  std::size_t produced_rows = 0;
  std::size_t failed_rows = 0;
  std::size_t produced_rows_a = 0;
  std::size_t produced_rows_b = 0;
  std::size_t failed_rows_a = 0;
  std::size_t failed_rows_b = 0;
};

struct SharedHybridClassCount {
  std::size_t invocation = 0;
  std::size_t tile_index = 0;
  std::size_t rows_seen_by_hybrid = 0;
  std::size_t class_none_count = 0;
  std::size_t class_bch_hard_decoded_count = 0;
  std::size_t class_clean_count = 0;
  std::size_t class_parity_only_count = 0;
  std::size_t class_one_main_count = 0;
  std::size_t class_one_main_plus_parity_count = 0;
  std::size_t class_two_main_count = 0;
  std::size_t class_suspicious_count = 0;
  std::size_t class_hard_fail_count = 0;
  std::size_t deferred_candidate_count = 0;
  std::size_t deferred_priority_0_count = 0;
  std::size_t deferred_priority_1_count = 0;
  std::size_t deferred_priority_2_count = 0;
  std::size_t deferred_priority_3_count = 0;
  std::size_t deferred_reclaimed_to_hard_finish_count = 0;
};

struct SharedRowMapEntry {
  std::size_t invocation = 0;
  std::size_t tile_index = 0;
  std::size_t merged_row = 0;
  std::size_t stream_id = 0;
  std::size_t source_local_row = 0;
  std::size_t source_global_row = 0;
  bool early_stop_hit = false;
  SharedHybridClass hybrid_class = SharedHybridClass::None;
  SharedRowFinalTag final_tag = SharedRowFinalTag::SoftDecode;
  bool scheduled_for_soft = false;
  bool scheduled_for_hard = false;
  bool produced_row = false;
};

struct Observability {
  float shared_quant_clip = 0.0f;
  StreamQuantizationStats stream_a_quantization;
  StreamQuantizationStats stream_b_quantization;
  SharedCoreAggregateStats shared_core;
  std::vector<SharedTileSample> tile_samples;
  std::vector<SharedHybridClassCount> hybrid_class_counts;
  std::vector<SharedRowMapEntry> row_map;
};

struct Result {
  PipelineResult stream_a;
  PipelineResult stream_b;
  Observability observability;
};

Result run_two_stream_shared(const Config& config);

}  // namespace newcode::two_stream_shared
