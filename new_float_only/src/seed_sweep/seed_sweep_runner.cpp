#include "new_float_only/seed_sweep_runner.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

#include "new_float_only/io/ensure_dir.hpp"
#include "new_float_only/sweep_task_runner.hpp"
#include "new_float_only/utils/now_stamp.hpp"

namespace new_float_only {
namespace {

std::string resolve_seed_mode(const SeedSweepConfig& config) {
  return config.bitgen_seeds.empty() ? "sequential" : "paired_list";
}

void write_trial_csv_header(std::ofstream& csv) {
  csv << "run_id,trial_index,bitgen_seed,channel_seed,"
         "pre_ber,pre_errors,pre_total,post_ber,post_errors,post_total\n";
}

void write_summary_csv_header(std::ofstream& csv) {
  csv << "run_id,label,seed_mode,ebn0_db,bits_per_symbol,normalize_extrinsic,"
         "generate_random_bits,max_workers,trials_requested,trials_completed,"
         "bitgen_seed_base,channel_seed_base,"
         "pre_ber,pre_errors,pre_total,pre_frame_error_trials,"
         "post_ber,post_errors,post_total,post_frame_error_trials\n";
}

void write_trial_csv_row(std::ofstream& csv,
                         const std::string& run_id,
                         const SeedTrialResult& trial) {
  csv << run_id << ","
      << trial.trial_index << ","
      << trial.bitgen_seed << ","
      << trial.channel_seed << ","
      << trial.pre_fec.ber << ","
      << trial.pre_fec.errors << ","
      << trial.pre_fec.total << ","
      << trial.post_fec.ber << ","
      << trial.post_fec.errors << ","
      << trial.post_fec.total << "\n";
}

void write_summary_csv_row(std::ofstream& csv,
                           const std::string& run_id,
                           const SeedSweepConfig& config,
                           const std::string& seed_mode,
                           unsigned max_workers,
                           const SeedSweepResult& result) {
  csv << run_id << ","
      << '"' << config.label << '"' << ","
      << seed_mode << ","
      << config.ebn0_db << ","
      << config.bits_per_symbol << ","
      << (config.normalize_extrinsic ? 1 : 0) << ","
      << (config.generate_random_bits ? 1 : 0) << ","
      << max_workers << ","
      << result.trials_requested << ","
      << result.trials_completed << ","
      << config.bitgen_seed_base << ","
      << config.channel_seed_base << ","
      << result.pre_fec.ber << ","
      << result.pre_fec.errors << ","
      << result.pre_fec.total << ","
      << result.pre_frame_error_trials << ","
      << result.post_fec.ber << ","
      << result.post_fec.errors << ","
      << result.post_fec.total << ","
      << result.post_frame_error_trials << "\n";
}

}  // namespace

SeedSweepResult run_seed_sweep(const SeedSweepConfig& config) {
  const Params params = normalize_sweep_decoder_config(config.decoder);
  const std::vector<SweepSeedPair> schedule =
      build_sweep_seed_schedule(config.trial_count,
                                config.bitgen_seed_base,
                                config.channel_seed_base,
                                config.bitgen_seeds,
                                config.channel_seeds);
  const std::string seed_mode = resolve_seed_mode(config);

  SweepPattern pattern;
  pattern.pattern_index = 0;
  pattern.pattern_label = config.label;
  pattern.alpha_list = params.ALPHA_LIST;
  pattern.beta_list = params.beta_list;
  const std::vector<SweepPattern> patterns = {pattern};
  const std::vector<SweepTask> tasks = build_sweep_tasks(patterns, schedule);

  SweepTaskRunnerConfig task_config;
  task_config.label = config.label;
  task_config.ebn0_db = config.ebn0_db;
  task_config.bits_per_symbol = config.bits_per_symbol;
  task_config.generate_random_bits = config.generate_random_bits;
  task_config.normalize_extrinsic = config.normalize_extrinsic;
  task_config.base_decoder = params;
  task_config.max_workers_override = config.max_workers_override;
  task_config.quiet_pipeline = config.quiet_pipeline;
  task_config.quiet_logs = config.quiet_logs;
  const unsigned max_workers = resolve_sweep_worker_count(task_config);

  if (!config.quiet_logs) {
    std::cout << "[SEED_SWEEP] label=" << config.label
              << ", seed_mode=" << seed_mode
              << ", trials=" << schedule.size()
              << ", workers=" << max_workers << "\n";
  }

  const std::vector<SweepTaskResult> task_results = run_sweep_tasks(task_config, tasks);
  const std::vector<SweepPatternSummary> summaries =
      aggregate_sweep_task_results(patterns, task_results);
  if (summaries.size() != 1) {
    throw std::runtime_error("Seed sweep expected exactly one pattern summary");
  }

  SeedSweepResult result;
  result.trials_requested = schedule.size();
  result.trials_completed = summaries.front().trials_completed;
  result.pre_fec = summaries.front().pre_fec;
  result.post_fec = summaries.front().post_fec;
  result.pre_frame_error_trials = summaries.front().pre_frame_error_trials;
  result.post_frame_error_trials = summaries.front().post_frame_error_trials;
  result.trial_results.reserve(task_results.size());
  for (const SweepTaskResult& task_result : task_results) {
    SeedTrialResult trial;
    trial.trial_index = task_result.seed_index;
    trial.bitgen_seed = task_result.bitgen_seed;
    trial.channel_seed = task_result.channel_seed;
    trial.pre_fec = task_result.pre_fec;
    trial.post_fec = task_result.post_fec;
    result.trial_results.push_back(trial);
  }

  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();

  std::ofstream summary_csv;
  if (config.write_summary_csv) {
    result.summary_csv_path =
        (data_dir / ("seed_sweep_results_" + run_id + ".csv")).string();
    summary_csv.open(result.summary_csv_path, std::ios::out | std::ios::trunc);
    if (!summary_csv.is_open()) {
      throw std::runtime_error("Failed to open summary CSV: " + result.summary_csv_path);
    }
    write_summary_csv_header(summary_csv);
    write_summary_csv_row(summary_csv, run_id, config, seed_mode, max_workers, result);
  }

  std::ofstream trial_csv;
  if (config.write_trial_csv) {
    result.trial_csv_path =
        (data_dir / ("seed_sweep_trials_" + run_id + ".csv")).string();
    trial_csv.open(result.trial_csv_path, std::ios::out | std::ios::trunc);
    if (!trial_csv.is_open()) {
      throw std::runtime_error("Failed to open trial CSV: " + result.trial_csv_path);
    }
    write_trial_csv_header(trial_csv);
    for (const SeedTrialResult& trial : result.trial_results) {
      write_trial_csv_row(trial_csv, run_id, trial);
    }
  }

  if (!config.quiet_logs) {
    std::cout << "[SEED_SWEEP] done: pre-BER=" << result.pre_fec.ber
              << " (" << result.pre_fec.errors << "/" << result.pre_fec.total << ")"
              << ", post-BER=" << result.post_fec.ber
              << " (" << result.post_fec.errors << "/" << result.post_fec.total << ")\n";
    std::cout << "[SEED_SWEEP] frame error trials: pre=" << result.pre_frame_error_trials
              << ", post=" << result.post_frame_error_trials << "\n";
    if (!result.summary_csv_path.empty()) {
      std::cout << "[SEED_SWEEP] summary CSV: " << result.summary_csv_path << "\n";
    }
    if (!result.trial_csv_path.empty()) {
      std::cout << "[SEED_SWEEP] trial CSV: " << result.trial_csv_path << "\n";
    }
  }

  return result;
}

}  // namespace new_float_only
