#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "new_float_only/params.hpp"
#include "new_float_only/rx/ber/ber.hpp"

namespace new_float_only {

/**
 * 单个 seed trial 的统计结果。
 * 每个 trial 固定对应一对 bitgen/channel seed，并返回该次链路的 pre/post BER。
 */
struct SeedTrialResult {
  std::size_t trial_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
  BerStats pre_fec;
  BerStats post_fec;
};

/**
 * 固定参数下的并行 Monte Carlo 配置。
 * 该配置不做参数网格扫描，只负责在相同 decoder/channel 参数下更换 seed 并聚合 BER。
 */
struct SeedSweepConfig {
  std::string label = "seed_sweep_float";
  float ebn0_db = 3.24f;
  unsigned bits_per_symbol = 2;
  bool generate_random_bits = true;
  bool normalize_extrinsic = true;
  DecoderConfig decoder{};

  std::size_t trial_count = 0;
  int bitgen_seed_base = 56456;
  int channel_seed_base = 57112;
  std::vector<int> bitgen_seeds;
  std::vector<int> channel_seeds;

  unsigned max_workers_override = 0;
  bool quiet_pipeline = true;
  bool quiet_logs = false;
  bool write_summary_csv = true;
  bool write_trial_csv = true;
};

/**
 * 整次 seed sweep 的最终结果。
 * trial_results 会按 trial_index 排序，便于调用方直接复核聚合前后的 BER。
 */
struct SeedSweepResult {
  std::size_t trials_requested = 0;
  std::size_t trials_completed = 0;
  BerStats pre_fec;
  BerStats post_fec;
  std::size_t pre_frame_error_trials = 0;
  std::size_t post_frame_error_trials = 0;
  std::vector<SeedTrialResult> trial_results;
  std::string summary_csv_path;
  std::string trial_csv_path;
};

/**
 * 并行运行固定参数的多 seed Monte Carlo 实验。
 * 输入：
 * 1. config 指定统一的链路参数、seed 调度规则、并行度和 CSV 输出选项；
 * 2. 当显式 seed 列表为空时，会使用 base+i 方式自动生成 seed；
 * 3. 当显式 seed 列表非空时，会按索引一一配对，不做全组合。
 * 输出：
 * 1. 返回聚合后的 pre/post BER；
 * 2. 同时回传每个 trial 的独立统计结果；
 * 3. 若开启 CSV 输出，会把路径写入 summary_csv_path/trial_csv_path。
 */
SeedSweepResult run_seed_sweep(const SeedSweepConfig& config);

}  // namespace new_float_only
