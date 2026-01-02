#include "tpc_sweep_detail.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>

#include "newcode/io/ensure_dir.hpp"
#include "newcode/utils/now_stamp.hpp"
namespace tpc_sweep {
namespace {

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

std::vector<int> resolve_seeds(const std::vector<int>& provided,
                               int count,
                               int fallback) {
  if (!provided.empty()) {
    return provided;
  }
  auto generated = detail::generate_random_seeds(count);
  if (generated.empty()) {
    generated.push_back(fallback);
  }
  return generated;
}

}  // namespace

namespace detail {

unsigned resolve_worker_count(const SweepParameterConfig& /*config*/) {
  unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  unsigned max_workers = std::max(1u, (hw * 3) / 4);
  if (const char* env = std::getenv("NTHREADS")) {
    try {
      int from_env = std::stoi(env);
      if (from_env > 0) {
        max_workers = static_cast<unsigned>(from_env);
      }
    } catch (...) {
      // ignore invalid environment variable
    }
  }
  return std::max(1u, max_workers);
}

Semaphore::Semaphore(std::size_t count)
  : count_(count ? count : 1) {}

void Semaphore::acquire() {
  std::unique_lock<std::mutex> lk(mutex_);
  cv_.wait(lk, [&] { return count_ > 0; });
  --count_;
}

void Semaphore::release() {
  std::lock_guard<std::mutex> lk(mutex_);
  ++count_;
  cv_.notify_one();
}

std::vector<ScenarioOutput> run_scenarios_parallel(
    const std::vector<SweepScenario>& scenarios,
    const SweepParameterConfig& config,
    unsigned max_workers_hint,
    const std::string& stage_tag,
    DualOut* log) {
  std::vector<ScenarioOutput> outputs;
  if (scenarios.empty()) {
    return outputs;
  }

  const unsigned max_workers =
      max_workers_hint > 0 ? max_workers_hint : resolve_worker_count(config);
  Semaphore sem(max_workers);

  newcode::Params base_params = config.base_params;
  base_params.BITGEN_RANDOM_BITS = config.generate_random_bits;
  newcode::tpc::TpcPipelineConfig pipeline_cfg;
  pipeline_cfg.bits_per_symbol = config.bits_per_symbol;
  pipeline_cfg.max_iters = config.max_iters;
  pipeline_cfg.num_blocks = config.num_blocks;
  pipeline_cfg.alpha_schedule = config.alpha_schedule;
  pipeline_cfg.beta_schedule = config.beta_schedule;
  pipeline_cfg.quiet = config.quiet_pipeline;

  const auto start_time = std::chrono::steady_clock::now();
  const std::size_t total_jobs = scenarios.size();
  std::atomic<std::size_t> finished{0};
  std::mutex progress_mutex;

  const std::string progress_prefix =
      stage_tag.empty() ? "[PROGRESS]" : ("[PROGRESS][" + stage_tag + "]");

  std::vector<std::future<ScenarioOutput>> futures;
  futures.reserve(scenarios.size());

  for (std::size_t idx = 0; idx < scenarios.size(); ++idx) {
    const auto& scenario = scenarios[idx];
    newcode::Params params = base_params;
    params.BITGEN_SEED = scenario.bitgen_seed;
    params.CHANNEL_SEED = scenario.channel_seed;

    sem.acquire();
    futures.emplace_back(std::async(
        std::launch::async,
        [idx,
         scenario,
         params,
         pipeline_cfg,
         &sem,
         log,
         start_time,
         &finished,
         total_jobs,
         &progress_mutex,
         progress_prefix,
         stage_tag]() -> ScenarioOutput {
          struct Releaser {
            detail::Semaphore& sem_ref;
            ~Releaser() { sem_ref.release(); }
          } releaser{sem};

          ScenarioOutput output;
          output.idx = idx;
          output.name = scenario.name;
          output.bitgen_seed = scenario.bitgen_seed;
          output.channel_seed = scenario.channel_seed;
          output.ebn0_db = scenario.ebn0_db;

          newcode::Params local_params = params;
          newcode::tpc::TpcPipelineConfig local_cfg = pipeline_cfg;
          const std::string label =
              stage_tag.empty() ? scenario.name : (scenario.name + "_" + stage_tag);
          output.result =
              newcode::tpc::run_tpc_pipeline(local_params, local_cfg, label, scenario.ebn0_db);
          output.ebn0_db = output.result.ebn0_db;

          const auto done = finished.fetch_add(1) + 1;
          const auto now = std::chrono::steady_clock::now();
          const auto elapsed = std::chrono::duration<double>(now - start_time);
          std::chrono::duration<double> eta(0.0);
          if (done < total_jobs && done > 0) {
            const double remaining_ratio =
                static_cast<double>(total_jobs - done) / static_cast<double>(done);
            eta = elapsed * remaining_ratio;
          }
          const double pct =
              static_cast<double>(done) * 100.0 / static_cast<double>(total_jobs);
          std::ostringstream oss;
          oss << progress_prefix << " " << done << "/" << total_jobs
              << " (" << std::fixed << std::setprecision(1) << pct << "%)"
              << " elapsed=" << format_duration(elapsed)
              << " ETA≈" << format_duration(eta);
          const std::string line = oss.str();
          std::lock_guard<std::mutex> lk(progress_mutex);
          if (log) {
            *log << line << "\n";
          } else {
            std::cout << line << "\n";
          }
          return output;
        }));
  }

  outputs.reserve(futures.size());
  for (auto& fut : futures) {
    try {
      outputs.emplace_back(fut.get());
    } catch (const std::exception& ex) {
      if (log) {
        *log << "[ERROR] worker threw: " << ex.what() << "\n";
      } else {
        std::cout << "[ERROR] worker threw: " << ex.what() << "\n";
      }
    } catch (...) {
      if (log) {
        *log << "[ERROR] worker threw unknown exception\n";
      } else {
        std::cout << "[ERROR] worker threw unknown exception\n";
      }
    }
  }

  if (outputs.empty()) {
    if (log) {
      *log << "[ERROR] No scenario completed successfully.\n";
    } else {
      std::cout << "[ERROR] No scenario completed successfully.\n";
    }
    return outputs;
  }

  std::sort(outputs.begin(), outputs.end(),
            [](const auto& a, const auto& b) { return a.idx < b.idx; });

  return outputs;
}

}  // namespace detail

int run_sweep(const SweepParameterConfig& config) {
  if (config.alpha_schedule.size() !=
      static_cast<std::size_t>(config.max_iters * 2)) {
    std::cout << "[ERROR] alpha_schedule length must be 2 * max_iters\n";
    return 1;
  }
  if (config.beta_schedule.size() !=
      static_cast<std::size_t>(config.max_iters * 2)) {
    std::cout << "[ERROR] beta_schedule length must be 2 * max_iters\n";
    return 1;
  }

  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();
  const std::string log_path = (data_dir / ("tpc_sweep_run_" + run_id + ".log")).string();
  const bool mirror_console = !config.quiet_logs;
  detail::DualOut out(std::cout, log_path, mirror_console);

  const std::string csv_path =
      (data_dir / ("tpc_sweep_results_" + run_id + ".csv")).string();
  detail::ensure_csv_header(csv_path);

  std::vector<int> bitgen_seeds =
      resolve_seeds(config.bitgen_seed_candidates,
                    config.bitgen_seed_count,
                    config.base_params.BITGEN_SEED);
  std::vector<int> channel_seeds =
      resolve_seeds(config.channel_seed_candidates,
                    config.channel_seed_count,
                    config.base_params.CHANNEL_SEED);
  const std::vector<float> ebn0_values = detail::build_ebn0_values(config);

  auto scenarios = detail::build_scenarios(config, ebn0_values, bitgen_seeds, channel_seeds);
  if (scenarios.empty()) {
    std::cout << "[ERROR] No scenarios generated.\n";
    return 1;
  }

  const unsigned max_workers = detail::resolve_worker_count(config);
  out << "[INFO] total scenarios = " << scenarios.size() << "\n";
  out << "[INFO] using up to " << max_workers << " workers\n";
  auto results = detail::run_scenarios_parallel(scenarios, config, max_workers, "", &out);
  if (results.empty()) {
    return 1;
  }

  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(8);

  for (const auto& pack : results) {
    const auto& result = pack.result;
    out << "[SUMMARY] " << pack.name
        << " Pre-FEC BER=" << result.pre_fec.ber
        << " (errs=" << result.pre_fec.errors << "/" << result.pre_fec.total << ")"
        << " | Post-FEC BER=" << result.post_fec.ber
        << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")"
        << " | Eb/N0=" << pack.ebn0_db
        << " | Seeds(bit/channel)=" << pack.bitgen_seed
        << "/" << pack.channel_seed << "\n";

    detail::write_csv_row(csv, utils::now_stamp(), run_id, scenarios[pack.idx], result, config);
  }
  return 0;
}

}  // namespace tpc_sweep
