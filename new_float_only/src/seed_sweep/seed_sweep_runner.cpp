#include "new_float_only/seed_sweep_runner.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#include "new_float_only/io/ensure_dir.hpp"
#include "new_float_only/pipeline_runner.hpp"
#include "new_float_only/utils/now_stamp.hpp"

namespace new_float_only {
namespace {

struct SeedPair {
  std::size_t trial_index = 0;
  int bitgen_seed = 0;
  int channel_seed = 0;
};

/**
 * 简单计数信号量。
 * sweep 会先 acquire 再启动 async worker，worker 结束时用析构自动 release，
 * 从而把并发中的试验数限制在 max_workers 以内。
 */
class Semaphore {
 public:
  explicit Semaphore(std::size_t count) : count_(count == 0 ? 1 : count) {}

  void acquire() {
    std::unique_lock<std::mutex> lock(mutex_);
    cv_.wait(lock, [&] { return count_ > 0; });
    --count_;
  }

  void release() {
    std::lock_guard<std::mutex> lock(mutex_);
    ++count_;
    cv_.notify_one();
  }

 private:
  std::mutex mutex_;
  std::condition_variable cv_;
  std::size_t count_;
};

std::string format_duration(std::chrono::duration<double> duration) {
  if (duration.count() < 0.0) {
    duration = std::chrono::duration<double>(0.0);
  }
  using Seconds = std::chrono::seconds;
  const auto total_seconds = std::chrono::duration_cast<Seconds>(duration).count();
  const auto hours = total_seconds / 3600;
  const auto minutes = (total_seconds % 3600) / 60;
  const auto seconds = total_seconds % 60;

  std::ostringstream oss;
  oss << std::setfill('0');
  if (hours > 0) {
    oss << hours << ':' << std::setw(2) << minutes << ':' << std::setw(2) << seconds;
  } else {
    oss << minutes << ':' << std::setw(2) << seconds;
  }
  return oss.str();
}

/**
 * 规范化解码参数，并在进入 sweep 前补齐 alpha/beta 列表。
 * 这样每个 worker 只需要复制已经稳定的 Params，不需要重复做参数修正。
 */
Params normalize_decoder_config(const SeedSweepConfig& config) {
  Params params = config.decoder;
  params.CHASE_NTEST = std::max(params.CHASE_NTEST, 1 << params.CHASE_L);

  if (!params.valid()) {
    throw std::invalid_argument("DecoderConfig is invalid");
  }
  if (!params.ALPHA_LIST.empty() && params.ALPHA_LIST.size() != params.TILES_PER_WIN) {
    throw std::invalid_argument("ALPHA_LIST size must equal TILES_PER_WIN");
  }
  if (!params.beta_list.empty() && params.beta_list.size() != params.TILES_PER_WIN) {
    throw std::invalid_argument("beta_list size must equal TILES_PER_WIN");
  }
  if (!params.ALPHA_LIST.empty()) {
    params.ALPHA = params.ALPHA_LIST.front();
  } else {
    params.ALPHA_LIST.assign(params.TILES_PER_WIN, params.ALPHA);
  }
  if (!params.beta_list.empty()) {
    params.beta = params.beta_list.front();
  } else {
    params.beta_list.assign(params.TILES_PER_WIN, params.beta);
  }
  return params;
}

void validate_sweep_side_effects(const Params& params) {
  if (params.DUMP_WORK_LLR) {
    throw std::invalid_argument(
        "Seed sweep does not support DUMP_WORK_LLR because parallel trials would overwrite the same file");
  }
  if (params.debug_trace.enable) {
    throw std::invalid_argument(
        "Seed sweep does not support debug_trace.enable because parallel trials would mix trace outputs");
  }
}

unsigned resolve_worker_count(const SeedSweepConfig& config) {
  unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  unsigned max_workers = std::max(1u, (hw * 3) / 4);
  if (config.max_workers_override > 0) {
    max_workers = config.max_workers_override;
  }
  if (const char* env = std::getenv("NTHREADS")) {
    try {
      const int from_env = std::stoi(env);
      if (from_env > 0) {
        max_workers = static_cast<unsigned>(from_env);
      }
    } catch (...) {
      // 忽略非法环境变量，回退到上面解析出的并行度。
    }
  }
  return std::max(1u, max_workers);
}

/**
 * 根据配置生成一一配对的 seed 调度表。
 * 两种模式：
 * 1. 显式列表模式：bitgen/channel 列表长度必须相等，trial_count=0 时自动取列表长度；
 * 2. 递增模式：使用 base+i 生成 seed，trial_count 必须大于 0。
 */
std::vector<SeedPair> build_seed_schedule(const SeedSweepConfig& config) {
  const bool has_bitgen_list = !config.bitgen_seeds.empty();
  const bool has_channel_list = !config.channel_seeds.empty();
  if (has_bitgen_list != has_channel_list) {
    throw std::invalid_argument(
        "bitgen_seeds and channel_seeds must either both be empty or both be provided");
  }

  std::vector<SeedPair> schedule;
  if (has_bitgen_list) {
    if (config.bitgen_seeds.size() != config.channel_seeds.size()) {
      throw std::invalid_argument("bitgen_seeds and channel_seeds must have the same length");
    }
    const std::size_t resolved_trials =
        (config.trial_count == 0) ? config.bitgen_seeds.size() : config.trial_count;
    if (resolved_trials != config.bitgen_seeds.size()) {
      throw std::invalid_argument(
          "trial_count must be 0 or equal to the explicit seed list length");
    }
    schedule.reserve(resolved_trials);
    for (std::size_t i = 0; i < resolved_trials; ++i) {
      schedule.push_back(SeedPair{
          .trial_index = i,
          .bitgen_seed = config.bitgen_seeds[i],
          .channel_seed = config.channel_seeds[i],
      });
    }
    return schedule;
  }

  if (config.trial_count == 0) {
    throw std::invalid_argument(
        "trial_count must be greater than 0 when explicit seed lists are not provided");
  }
  schedule.reserve(config.trial_count);
  for (std::size_t i = 0; i < config.trial_count; ++i) {
    schedule.push_back(SeedPair{
        .trial_index = i,
        .bitgen_seed = config.bitgen_seed_base + static_cast<int>(i),
        .channel_seed = config.channel_seed_base + static_cast<int>(i),
    });
  }
  return schedule;
}

PipelineConfig build_pipeline_config(const SeedSweepConfig& config,
                                     int bitgen_seed,
                                     int channel_seed) {
  PipelineConfig pipeline;
  pipeline.normalize_extrinsic = config.normalize_extrinsic;
  pipeline.bits_per_symbol = config.bits_per_symbol;
  pipeline.bitgen_seed = bitgen_seed;
  pipeline.channel_seed = channel_seed;
  pipeline.generate_random_bits = config.generate_random_bits;
  pipeline.collect_error_positions = false;
  pipeline.quiet = config.quiet_pipeline;
  return pipeline;
}

BerStats aggregate_stats(std::size_t errors, std::size_t total) {
  BerStats stats;
  stats.errors = errors;
  stats.total = total;
  stats.ber = (total == 0) ? 0.0 : static_cast<double>(errors) / static_cast<double>(total);
  return stats;
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

void print_overview(const SeedSweepConfig& config,
                    const Params& params,
                    const std::string& seed_mode,
                    std::size_t trial_count,
                    unsigned max_workers) {
  if (config.quiet_logs) {
    return;
  }
  std::cout << "[SEED_SWEEP] label=" << config.label
            << ", Eb/N0=" << config.ebn0_db
            << " dB, trials=" << trial_count
            << ", workers=" << max_workers << "\n";
  std::cout << "[SEED_SWEEP] bits_per_symbol=" << config.bits_per_symbol
            << ", normalize_extrinsic=" << (config.normalize_extrinsic ? "ON" : "OFF")
            << ", generate_random_bits=" << (config.generate_random_bits ? "ON" : "OFF")
            << ", seed_mode=" << seed_mode << "\n";
  std::cout << "[SEED_SWEEP] CHASE_L=" << params.CHASE_L
            << ", CHASE_NTEST=" << params.CHASE_NTEST
            << ", TILES_PER_WIN=" << params.TILES_PER_WIN
            << ", quiet_pipeline=" << (config.quiet_pipeline ? "ON" : "OFF") << "\n";
}

void print_progress(const SeedSweepConfig& config,
                    std::size_t done,
                    std::size_t total,
                    std::chrono::steady_clock::time_point start_time,
                    std::mutex& progress_mutex) {
  if (config.quiet_logs) {
    return;
  }

  const auto now = std::chrono::steady_clock::now();
  const auto elapsed = std::chrono::duration<double>(now - start_time);
  std::chrono::duration<double> eta(0.0);
  if (done < total && done > 0) {
    const double remaining_ratio =
        static_cast<double>(total - done) / static_cast<double>(done);
    eta = elapsed * remaining_ratio;
  }
  const double pct = static_cast<double>(done) * 100.0 / static_cast<double>(total);

  std::ostringstream oss;
  oss << "[PROGRESS] " << done << "/" << total
      << " (" << std::fixed << std::setprecision(1) << pct << "%)"
      << " elapsed=" << format_duration(elapsed)
      << " ETA≈" << format_duration(eta);

  std::lock_guard<std::mutex> lock(progress_mutex);
  std::cout << oss.str() << "\n";
}

}  // namespace

SeedSweepResult run_seed_sweep(const SeedSweepConfig& config) {
  const Params params = normalize_decoder_config(config);
  validate_sweep_side_effects(params);

  const std::vector<SeedPair> schedule = build_seed_schedule(config);
  const std::string seed_mode =
      schedule.empty() ? "empty" : (!config.bitgen_seeds.empty() ? "paired_list" : "sequential");
  const unsigned max_workers = resolve_worker_count(config);

  print_overview(config, params, seed_mode, schedule.size(), max_workers);

  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();

  SeedSweepResult result;
  result.trials_requested = schedule.size();

  std::ofstream summary_csv;
  if (config.write_summary_csv) {
    result.summary_csv_path =
        (data_dir / ("seed_sweep_results_" + run_id + ".csv")).string();
    summary_csv.open(result.summary_csv_path, std::ios::out | std::ios::trunc);
    if (!summary_csv.is_open()) {
      throw std::runtime_error("Failed to open summary CSV: " + result.summary_csv_path);
    }
    write_summary_csv_header(summary_csv);
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
  }

  Semaphore sem(max_workers);
  std::atomic<std::size_t> finished{0};
  std::mutex progress_mutex;
  const auto start_time = std::chrono::steady_clock::now();

  std::vector<std::future<SeedTrialResult>> futures;
  futures.reserve(schedule.size());
  const std::size_t total_trials = schedule.size();
  for (const SeedPair& pair : schedule) {
    sem.acquire();
    futures.emplace_back(std::async(std::launch::async,
                                    [pair,
                                     params,
                                     &config,
                                     &sem,
                                     &finished,
                                     &progress_mutex,
                                     start_time,
                                     total_trials]() -> SeedTrialResult {
      struct Releaser {
        Semaphore& sem_ref;
        ~Releaser() { sem_ref.release(); }
      } releaser{sem};

      const PipelineConfig pipeline_config =
          build_pipeline_config(config, pair.bitgen_seed, pair.channel_seed);
      const std::string label = config.label + "_trial" + std::to_string(pair.trial_index);
      const PipelineResult pipeline_result =
          run_pipeline(params, pipeline_config, label, config.ebn0_db);

      SeedTrialResult trial;
      trial.trial_index = pair.trial_index;
      trial.bitgen_seed = pair.bitgen_seed;
      trial.channel_seed = pair.channel_seed;
      trial.pre_fec = pipeline_result.pre_fec;
      trial.post_fec = pipeline_result.post_fec;

      const std::size_t done = finished.fetch_add(1) + 1;
      print_progress(config, done, total_trials, start_time, progress_mutex);
      return trial;
    }));
  }

  result.trial_results.reserve(futures.size());
  std::exception_ptr worker_error;
  for (auto& future : futures) {
    try {
      result.trial_results.push_back(future.get());
    } catch (...) {
      if (!worker_error) {
        worker_error = std::current_exception();
      }
    }
  }
  if (worker_error) {
    std::rethrow_exception(worker_error);
  }

  std::sort(result.trial_results.begin(),
            result.trial_results.end(),
            [](const SeedTrialResult& lhs, const SeedTrialResult& rhs) {
              return lhs.trial_index < rhs.trial_index;
            });

  std::size_t pre_errors = 0;
  std::size_t pre_total = 0;
  std::size_t post_errors = 0;
  std::size_t post_total = 0;
  for (const SeedTrialResult& trial : result.trial_results) {
    pre_errors += trial.pre_fec.errors;
    pre_total += trial.pre_fec.total;
    post_errors += trial.post_fec.errors;
    post_total += trial.post_fec.total;
    result.pre_frame_error_trials += (trial.pre_fec.errors > 0) ? 1u : 0u;
    result.post_frame_error_trials += (trial.post_fec.errors > 0) ? 1u : 0u;
    if (trial_csv.is_open()) {
      write_trial_csv_row(trial_csv, run_id, trial);
    }
  }

  result.trials_completed = result.trial_results.size();
  result.pre_fec = aggregate_stats(pre_errors, pre_total);
  result.post_fec = aggregate_stats(post_errors, post_total);

  if (summary_csv.is_open()) {
    write_summary_csv_row(summary_csv, run_id, config, seed_mode, max_workers, result);
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
