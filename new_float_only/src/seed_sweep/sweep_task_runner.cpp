#include "new_float_only/sweep_task_runner.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <future>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "new_float_only/pipeline_runner.hpp"

namespace new_float_only {
namespace {

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

void validate_parallel_sweep_side_effects(const Params& params) {
  if (params.DUMP_WORK_LLR) {
    throw std::invalid_argument(
        "Parallel sweep does not support DUMP_WORK_LLR because tasks would overwrite the same file");
  }
  if (params.debug_trace.enable) {
    throw std::invalid_argument(
        "Parallel sweep does not support debug_trace.enable because tasks would mix trace outputs");
  }
}

unsigned resolve_worker_count_impl(const SweepTaskRunnerConfig& config) {
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
    }
  }
  return std::max(1u, max_workers);
}

PipelineConfig build_pipeline_config(const SweepTaskRunnerConfig& config,
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

std::string resolve_pattern_label(const SweepPattern& pattern) {
  if (!pattern.pattern_label.empty()) {
    return pattern.pattern_label;
  }
  return "pattern_" + std::to_string(pattern.pattern_index);
}

void validate_task_lists(const Params& params, const SweepTask& task) {
  if (!task.alpha_list.empty() && task.alpha_list.size() != params.TILES_PER_WIN) {
    throw std::invalid_argument("Task alpha_list size must equal TILES_PER_WIN");
  }
  if (!task.beta_list.empty() && task.beta_list.size() != params.TILES_PER_WIN) {
    throw std::invalid_argument("Task beta_list size must equal TILES_PER_WIN");
  }
}

Params build_task_params(const Params& base_params, const SweepTask& task) {
  Params params = base_params;
  validate_task_lists(params, task);

  if (!task.alpha_list.empty()) {
    params.ALPHA_LIST = task.alpha_list;
    params.ALPHA = task.alpha_list.front();
  }
  if (!task.beta_list.empty()) {
    params.beta_list = task.beta_list;
    params.beta = task.beta_list.front();
  }
  return params;
}

void print_overview(const SweepTaskRunnerConfig& config,
                    std::size_t task_count,
                    unsigned max_workers) {
  if (config.quiet_logs) {
    return;
  }

  std::cout << "[SWEEP_TASK] label=" << config.label
            << ", Eb/N0=" << config.ebn0_db
            << " dB, tasks=" << task_count
            << ", workers=" << max_workers << "\n";
  std::cout << "[SWEEP_TASK] bits_per_symbol=" << config.bits_per_symbol
            << ", normalize_extrinsic=" << (config.normalize_extrinsic ? "ON" : "OFF")
            << ", generate_random_bits=" << (config.generate_random_bits ? "ON" : "OFF")
            << ", quiet_pipeline=" << (config.quiet_pipeline ? "ON" : "OFF") << "\n";
}

void print_progress(const SweepTaskRunnerConfig& config,
                    std::size_t done,
                    std::size_t total,
                    std::chrono::steady_clock::time_point start_time,
                    std::mutex& progress_mutex) {
  if (config.quiet_logs || total == 0) {
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

Params normalize_sweep_decoder_config(const DecoderConfig& decoder) {
  Params params = decoder;
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

unsigned resolve_sweep_worker_count(const SweepTaskRunnerConfig& config) {
  return resolve_worker_count_impl(config);
}

std::vector<SweepSeedPair> build_sweep_seed_schedule(std::size_t trial_count,
                                                     int bitgen_seed_base,
                                                     int channel_seed_base,
                                                     const std::vector<int>& bitgen_seeds,
                                                     const std::vector<int>& channel_seeds) {
  const bool has_bitgen_list = !bitgen_seeds.empty();
  const bool has_channel_list = !channel_seeds.empty();
  if (has_bitgen_list != has_channel_list) {
    throw std::invalid_argument(
        "bitgen_seeds and channel_seeds must either both be empty or both be provided");
  }

  std::vector<SweepSeedPair> schedule;
  if (has_bitgen_list) {
    if (bitgen_seeds.size() != channel_seeds.size()) {
      throw std::invalid_argument("bitgen_seeds and channel_seeds must have the same length");
    }
    const std::size_t resolved_trials = (trial_count == 0) ? bitgen_seeds.size() : trial_count;
    if (resolved_trials != bitgen_seeds.size()) {
      throw std::invalid_argument(
          "trial_count must be 0 or equal to the explicit seed list length");
    }

    schedule.reserve(resolved_trials);
    for (std::size_t i = 0; i < resolved_trials; ++i) {
      schedule.push_back(SweepSeedPair{
          .seed_index = i,
          .bitgen_seed = bitgen_seeds[i],
          .channel_seed = channel_seeds[i],
      });
    }
    return schedule;
  }

  if (trial_count == 0) {
    throw std::invalid_argument(
        "trial_count must be greater than 0 when explicit seed lists are not provided");
  }

  schedule.reserve(trial_count);
  for (std::size_t i = 0; i < trial_count; ++i) {
    schedule.push_back(SweepSeedPair{
        .seed_index = i,
        .bitgen_seed = bitgen_seed_base + static_cast<int>(i),
        .channel_seed = channel_seed_base + static_cast<int>(i),
    });
  }
  return schedule;
}

std::vector<SweepTask> build_sweep_tasks(const std::vector<SweepPattern>& patterns,
                                         const std::vector<SweepSeedPair>& seed_schedule) {
  std::vector<SweepTask> tasks;
  tasks.reserve(patterns.size() * seed_schedule.size());

  std::size_t task_index = 0;
  for (const SweepPattern& pattern : patterns) {
    const std::string pattern_label = resolve_pattern_label(pattern);
    for (const SweepSeedPair& seed_pair : seed_schedule) {
      tasks.push_back(SweepTask{
          .task_index = task_index++,
          .pattern_index = pattern.pattern_index,
          .seed_index = seed_pair.seed_index,
          .pattern_label = pattern_label,
          .alpha_list = pattern.alpha_list,
          .beta_list = pattern.beta_list,
          .bitgen_seed = seed_pair.bitgen_seed,
          .channel_seed = seed_pair.channel_seed,
      });
    }
  }

  return tasks;
}

std::vector<SweepTaskResult> run_sweep_tasks(const SweepTaskRunnerConfig& config,
                                             const std::vector<SweepTask>& tasks) {
  validate_parallel_sweep_side_effects(config.base_decoder);
  const unsigned max_workers = resolve_worker_count_impl(config);
  print_overview(config, tasks.size(), max_workers);

  Semaphore sem(max_workers);
  std::atomic<std::size_t> finished{0};
  std::mutex progress_mutex;
  const auto start_time = std::chrono::steady_clock::now();

  std::vector<std::future<SweepTaskResult>> futures;
  futures.reserve(tasks.size());
  const std::size_t total_tasks = tasks.size();
  for (const SweepTask& task : tasks) {
    sem.acquire();
    futures.emplace_back(std::async(std::launch::async,
                                    [task,
                                     &config,
                                     &sem,
                                     &finished,
                                     &progress_mutex,
                                     start_time,
                                     total_tasks]() -> SweepTaskResult {
      struct Releaser {
        Semaphore& sem_ref;
        ~Releaser() { sem_ref.release(); }
      } releaser{sem};

      const Params task_params = build_task_params(config.base_decoder, task);
      const PipelineConfig pipeline_config =
          build_pipeline_config(config, task.bitgen_seed, task.channel_seed);
      const std::string task_label =
          config.label + "_p" + std::to_string(task.pattern_index) +
          "_s" + std::to_string(task.seed_index);
      const PipelineResult pipeline_result =
          run_pipeline(task_params, pipeline_config, task_label, config.ebn0_db);

      SweepTaskResult result;
      result.task_index = task.task_index;
      result.pattern_index = task.pattern_index;
      result.seed_index = task.seed_index;
      result.pattern_label = task.pattern_label;
      result.alpha_list = task.alpha_list;
      result.beta_list = task.beta_list;
      result.bitgen_seed = task.bitgen_seed;
      result.channel_seed = task.channel_seed;
      result.pre_fec = pipeline_result.pre_fec;
      result.post_fec = pipeline_result.post_fec;

      const std::size_t done = finished.fetch_add(1) + 1;
      print_progress(config, done, total_tasks, start_time, progress_mutex);
      return result;
    }));
  }

  std::vector<SweepTaskResult> results;
  results.reserve(futures.size());
  std::exception_ptr worker_error;
  for (auto& future : futures) {
    try {
      results.push_back(future.get());
    } catch (...) {
      if (!worker_error) {
        worker_error = std::current_exception();
      }
    }
  }
  if (worker_error) {
    std::rethrow_exception(worker_error);
  }

  std::sort(results.begin(),
            results.end(),
            [](const SweepTaskResult& lhs, const SweepTaskResult& rhs) {
              return lhs.task_index < rhs.task_index;
            });

  if (!config.quiet_logs) {
    std::cout << "[SWEEP_TASK] done: completed=" << results.size()
              << "/" << total_tasks << "\n";
  }
  return results;
}

std::vector<SweepPatternSummary> aggregate_sweep_task_results(
    const std::vector<SweepPattern>& patterns,
    const std::vector<SweepTaskResult>& task_results) {
  std::vector<SweepPatternSummary> summaries(patterns.size());
  for (const SweepPattern& pattern : patterns) {
    if (pattern.pattern_index >= summaries.size()) {
      throw std::invalid_argument("pattern_index must be a dense range starting from 0");
    }

    SweepPatternSummary& summary = summaries[pattern.pattern_index];
    summary.pattern_index = pattern.pattern_index;
    summary.pattern_label = resolve_pattern_label(pattern);
    summary.alpha_list = pattern.alpha_list;
    summary.beta_list = pattern.beta_list;
  }

  std::vector<std::size_t> pre_errors(summaries.size(), 0);
  std::vector<std::size_t> pre_total(summaries.size(), 0);
  std::vector<std::size_t> post_errors(summaries.size(), 0);
  std::vector<std::size_t> post_total(summaries.size(), 0);

  for (const SweepTaskResult& task_result : task_results) {
    if (task_result.pattern_index >= summaries.size()) {
      throw std::invalid_argument("task_result.pattern_index is out of range");
    }

    SweepPatternSummary& summary = summaries[task_result.pattern_index];
    summary.trials_completed += 1;
    summary.pre_frame_error_trials += (task_result.pre_fec.errors > 0) ? 1u : 0u;
    summary.post_frame_error_trials += (task_result.post_fec.errors > 0) ? 1u : 0u;
    pre_errors[task_result.pattern_index] += task_result.pre_fec.errors;
    pre_total[task_result.pattern_index] += task_result.pre_fec.total;
    post_errors[task_result.pattern_index] += task_result.post_fec.errors;
    post_total[task_result.pattern_index] += task_result.post_fec.total;
  }

  for (std::size_t i = 0; i < summaries.size(); ++i) {
    summaries[i].pre_fec = aggregate_stats(pre_errors[i], pre_total[i]);
    summaries[i].post_fec = aggregate_stats(post_errors[i], post_total[i]);
  }

  return summaries;
}

}  // namespace new_float_only
