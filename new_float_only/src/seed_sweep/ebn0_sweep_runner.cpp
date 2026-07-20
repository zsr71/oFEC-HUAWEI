#include "new_float_only/ebn0_sweep_runner.hpp"

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
#include "new_float_only/sweep_task_runner.hpp"
#include "new_float_only/utils/now_stamp.hpp"

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

struct Ebn0Task {
  std::size_t task_index = 0;
  std::size_t ebn0_index = 0;
  std::size_t trial_index = 0;
  float ebn0_db = 0.0f;
  int bitgen_seed = 0;
  int channel_seed = 0;
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

void validate_parallel_side_effects(const Params& params) {
  if (params.DUMP_WORK_LLR) {
    throw std::invalid_argument(
        "Parallel Eb/N0 sweep does not support DUMP_WORK_LLR because tasks would overwrite the same file");
  }
  if (params.debug_trace.enable) {
    throw std::invalid_argument(
        "Parallel Eb/N0 sweep does not support debug_trace.enable because tasks would mix trace outputs");
  }
}

unsigned resolve_worker_count(const Ebn0SweepConfig& config) {
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

PipelineConfig build_pipeline_config(const Ebn0SweepConfig& config,
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

std::vector<float> build_ebn0_values(float start, float end, int points) {
  if (points <= 0) {
    throw std::invalid_argument("ebn0_points must be greater than 0");
  }

  std::vector<float> values(static_cast<std::size_t>(points), start);
  if (points == 1) {
    return values;
  }

  for (int i = 0; i < points; ++i) {
    const float t = static_cast<float>(i) / static_cast<float>(points - 1);
    values[static_cast<std::size_t>(i)] = start + (end - start) * t;
  }
  return values;
}

std::vector<Ebn0Task> build_ebn0_tasks(const std::vector<float>& ebn0_values,
                                       const std::vector<SweepSeedPair>& seed_schedule) {
  std::vector<Ebn0Task> tasks;
  tasks.reserve(ebn0_values.size() * seed_schedule.size());

  std::size_t task_index = 0;
  for (std::size_t ebn0_index = 0; ebn0_index < ebn0_values.size(); ++ebn0_index) {
    for (const SweepSeedPair& seed_pair : seed_schedule) {
      tasks.push_back(Ebn0Task{
          .task_index = task_index++,
          .ebn0_index = ebn0_index,
          .trial_index = seed_pair.seed_index,
          .ebn0_db = ebn0_values[ebn0_index],
          .bitgen_seed = seed_pair.bitgen_seed,
          .channel_seed = seed_pair.channel_seed,
      });
    }
  }
  return tasks;
}

void print_overview(const Ebn0SweepConfig& config,
                    std::size_t point_count,
                    std::size_t task_count,
                    unsigned max_workers) {
  if (config.quiet_logs) {
    return;
  }

  std::cout << "[EBN0_SWEEP] label=" << config.label
            << ", Eb/N0 start=" << config.ebn0_start
            << ", end=" << config.ebn0_end
            << ", points=" << point_count
            << ", tasks=" << task_count
            << ", workers=" << max_workers << "\n";
  std::cout << "[EBN0_SWEEP] bits_per_symbol=" << config.bits_per_symbol
            << ", normalize_extrinsic=" << (config.normalize_extrinsic ? "ON" : "OFF")
            << ", generate_random_bits=" << (config.generate_random_bits ? "ON" : "OFF")
            << ", quiet_pipeline=" << (config.quiet_pipeline ? "ON" : "OFF") << "\n";
}

void print_progress(const Ebn0SweepConfig& config,
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

std::string format_float_list(const std::vector<float>& values) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ';';
    }
    oss << std::fixed << std::setprecision(3) << values[i];
  }
  return oss.str();
}

std::string format_int_list(const std::vector<int>& values) {
  std::ostringstream oss;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (i) {
      oss << ';';
    }
    oss << values[i];
  }
  return oss.str();
}

void write_summary_csv_header(std::ofstream& csv) {
  csv << "run_id,label,ebn0_start,ebn0_end,ebn0_points,ebn0_index,ebn0_db,"
         "bits_per_symbol,generate_random_bits,normalize_extrinsic,max_workers_override,"
         "trial_count,bitgen_seed_base,channel_seed_base,"
         "normalize_known_prefix_tail,chase_l,chase_ntest,hard_decode_default,hard_tile_list,"
         "pre_ber,pre_errors,pre_total,pre_frame_error_trials,"
         "post_ber,post_errors,post_total,post_frame_error_trials,"
         "alpha_list,beta_list\n";
}

void write_trial_csv_header(std::ofstream& csv) {
  csv << "run_id,label,ebn0_index,ebn0_db,trial_index,bitgen_seed,channel_seed,"
         "pre_ber,pre_errors,pre_total,post_ber,post_errors,post_total\n";
}

void write_summary_csv_row(std::ofstream& csv,
                           const std::string& run_id,
                           const Ebn0SweepConfig& config,
                           const Params& params,
                           const Ebn0SweepPointResult& point) {
  csv << run_id << ","
      << '"' << config.label << '"' << ","
      << config.ebn0_start << ","
      << config.ebn0_end << ","
      << config.ebn0_points << ","
      << point.ebn0_index << ","
      << point.ebn0_db << ","
      << config.bits_per_symbol << ","
      << (config.generate_random_bits ? 1 : 0) << ","
      << (config.normalize_extrinsic ? 1 : 0) << ","
      << config.max_workers_override << ","
      << config.trial_count << ","
      << config.bitgen_seed_base << ","
      << config.channel_seed_base << ","
      << (params.NORMALIZE_KNOWN_PREFIX_TAIL ? 1 : 0) << ","
      << params.CHASE_L << ","
      << params.CHASE_NTEST << ","
      << (params.HARD_DECODE_DEFAULT ? 1 : 0) << ","
      << '"' << format_int_list(params.HARD_TILE_LIST) << '"' << ","
      << point.trials_requested << ","
      << point.trials_completed << ","
      << point.pre_fec.ber << ","
      << point.pre_fec.errors << ","
      << point.pre_fec.total << ","
      << point.pre_frame_error_trials << ","
      << point.post_fec.ber << ","
      << point.post_fec.errors << ","
      << point.post_fec.total << ","
      << point.post_frame_error_trials << ","
      << '"' << format_float_list(params.ALPHA_LIST) << '"' << ","
      << '"' << format_float_list(params.beta_list) << '"' << "\n";
}

void write_trial_csv_row(std::ofstream& csv,
                         const std::string& run_id,
                         const Ebn0SweepConfig& config,
                         const Ebn0SweepTrialResult& trial) {
  csv << run_id << ","
      << '"' << config.label << '"' << ","
      << trial.ebn0_index << ","
      << trial.ebn0_db << ","
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

}  // namespace

Ebn0SweepResult run_ebn0_sweep(const Ebn0SweepConfig& config) {
  const Params params = normalize_sweep_decoder_config(config.decoder);
  validate_parallel_side_effects(params);

  const std::vector<float> ebn0_values =
      build_ebn0_values(config.ebn0_start, config.ebn0_end, config.ebn0_points);
  const std::vector<SweepSeedPair> seed_schedule =
      build_sweep_seed_schedule(config.trial_count,
                                config.bitgen_seed_base,
                                config.channel_seed_base,
                                config.bitgen_seeds,
                                config.channel_seeds);
  const std::vector<Ebn0Task> tasks = build_ebn0_tasks(ebn0_values, seed_schedule);
  const unsigned max_workers = resolve_worker_count(config);

  print_overview(config, ebn0_values.size(), tasks.size(), max_workers);

  Semaphore sem(max_workers);
  std::atomic<std::size_t> finished{0};
  std::mutex progress_mutex;
  const auto start_time = std::chrono::steady_clock::now();

  std::vector<std::future<Ebn0SweepTrialResult>> futures;
  futures.reserve(tasks.size());
  for (const Ebn0Task& task : tasks) {
    sem.acquire();
    futures.emplace_back(std::async(std::launch::async,
                                    [task,
                                     &config,
                                     &params,
                                     &sem,
                                     &finished,
                                     &progress_mutex,
                                     start_time,
                                     total_tasks = tasks.size()]() -> Ebn0SweepTrialResult {
      struct Releaser {
        Semaphore& sem_ref;
        ~Releaser() { sem_ref.release(); }
      } releaser{sem};

      const PipelineConfig pipeline_config =
          build_pipeline_config(config, task.bitgen_seed, task.channel_seed);
      const std::string task_label =
          config.label + "_e" + std::to_string(task.ebn0_index) +
          "_s" + std::to_string(task.trial_index);
      const PipelineResult pipeline_result =
          run_pipeline(params, pipeline_config, task_label, task.ebn0_db);

      Ebn0SweepTrialResult result;
      result.task_index = task.task_index;
      result.ebn0_index = task.ebn0_index;
      result.trial_index = task.trial_index;
      result.ebn0_db = task.ebn0_db;
      result.bitgen_seed = task.bitgen_seed;
      result.channel_seed = task.channel_seed;
      result.pre_fec = pipeline_result.pre_fec;
      result.post_fec = pipeline_result.post_fec;

      const std::size_t done = finished.fetch_add(1) + 1;
      print_progress(config, done, total_tasks, start_time, progress_mutex);
      return result;
    }));
  }

  Ebn0SweepResult result;
  result.points_requested = ebn0_values.size();
  result.ebn0_values = ebn0_values;
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
            [](const Ebn0SweepTrialResult& lhs, const Ebn0SweepTrialResult& rhs) {
              return lhs.task_index < rhs.task_index;
            });

  result.point_results.resize(ebn0_values.size());
  std::vector<std::size_t> pre_errors(ebn0_values.size(), 0);
  std::vector<std::size_t> pre_total(ebn0_values.size(), 0);
  std::vector<std::size_t> post_errors(ebn0_values.size(), 0);
  std::vector<std::size_t> post_total(ebn0_values.size(), 0);

  for (std::size_t i = 0; i < ebn0_values.size(); ++i) {
    result.point_results[i].ebn0_index = i;
    result.point_results[i].ebn0_db = ebn0_values[i];
    result.point_results[i].trials_requested = seed_schedule.size();
  }

  for (const Ebn0SweepTrialResult& trial : result.trial_results) {
    Ebn0SweepPointResult& point = result.point_results[trial.ebn0_index];
    point.trials_completed += 1;
    point.pre_frame_error_trials += (trial.pre_fec.errors > 0) ? 1u : 0u;
    point.post_frame_error_trials += (trial.post_fec.errors > 0) ? 1u : 0u;
    pre_errors[trial.ebn0_index] += trial.pre_fec.errors;
    pre_total[trial.ebn0_index] += trial.pre_fec.total;
    post_errors[trial.ebn0_index] += trial.post_fec.errors;
    post_total[trial.ebn0_index] += trial.post_fec.total;
  }

  for (std::size_t i = 0; i < result.point_results.size(); ++i) {
    result.point_results[i].pre_fec = aggregate_stats(pre_errors[i], pre_total[i]);
    result.point_results[i].post_fec = aggregate_stats(post_errors[i], post_total[i]);
    if (result.point_results[i].trials_completed > 0) {
      result.points_completed += 1;
    }
  }

  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();

  if (config.write_summary_csv) {
    result.summary_csv_path =
        (data_dir / ("ebn0_sweep_results_" + run_id + ".csv")).string();
    std::ofstream summary_csv(result.summary_csv_path, std::ios::out | std::ios::trunc);
    if (!summary_csv.is_open()) {
      throw std::runtime_error("Failed to open summary CSV: " + result.summary_csv_path);
    }
    write_summary_csv_header(summary_csv);
    for (const Ebn0SweepPointResult& point : result.point_results) {
      write_summary_csv_row(summary_csv, run_id, config, params, point);
    }
  }

  if (config.write_trial_csv) {
    result.trial_csv_path =
        (data_dir / ("ebn0_sweep_trials_" + run_id + ".csv")).string();
    std::ofstream trial_csv(result.trial_csv_path, std::ios::out | std::ios::trunc);
    if (!trial_csv.is_open()) {
      throw std::runtime_error("Failed to open trial CSV: " + result.trial_csv_path);
    }
    write_trial_csv_header(trial_csv);
    for (const Ebn0SweepTrialResult& trial : result.trial_results) {
      write_trial_csv_row(trial_csv, run_id, config, trial);
    }
  }

  if (!config.quiet_logs) {
    for (const Ebn0SweepPointResult& point : result.point_results) {
      std::cout << "[EBN0_SWEEP] Eb/N0=" << std::fixed << std::setprecision(3)
                << point.ebn0_db << " dB"
                << ", post-BER=" << point.post_fec.ber
                << " (" << point.post_fec.errors << "/" << point.post_fec.total << ")"
                << ", trials=" << point.trials_completed
                << "/" << point.trials_requested << "\n";
    }
    if (!result.summary_csv_path.empty()) {
      std::cout << "[EBN0_SWEEP] summary CSV: " << result.summary_csv_path << "\n";
    }
    if (!result.trial_csv_path.empty()) {
      std::cout << "[EBN0_SWEEP] trial CSV: " << result.trial_csv_path << "\n";
    }
  }

  return result;
}

}  // namespace new_float_only
