#include "ofec_sweep_detail.hpp"
#include "newcode/io/ensure_dir.hpp"
#include "newcode/ofec/mux/mux_config_validate.hpp"
#include "newcode/ofec/mux/mux_group_config_validate.hpp"
#include "newcode/utils/now_stamp.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <thread>

namespace ofec_sweep {
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

template <typename T>
std::string join_compact(const std::vector<T>& values) {
  std::ostringstream oss;
  oss << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    oss << values[i];
    if (i + 1 < values.size()) {
      oss << ",";
    }
  }
  oss << "]";
  return oss.str();
}

std::string join_compact(const std::vector<bool>& values) {
  std::ostringstream oss;
  oss << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    oss << (values[i] ? "true" : "false");
    if (i + 1 < values.size()) {
      oss << ",";
    }
  }
  oss << "]";
  return oss.str();
}

std::string join_compact(const std::vector<std::string>& values) {
  std::ostringstream oss;
  oss << "[";
  for (std::size_t i = 0; i < values.size(); ++i) {
    oss << values[i];
    if (i + 1 < values.size()) {
      oss << ",";
    }
  }
  oss << "]";
  return oss.str();
}

std::string summarize_patterns(
    const std::vector<ofec_sweep::ExplicitAlphaBetaPattern>& patterns) {
  std::ostringstream oss;
  oss << "[";
  for (std::size_t i = 0; i < patterns.size(); ++i) {
    const auto& p = patterns[i];
    oss << "{label=" << (p.label.empty() ? "<empty>" : p.label)
        << ",alpha=" << join_compact(p.alpha_list)
        << ",beta=" << join_compact(p.beta_list)
        << ",es_beta=" << join_compact(p.early_stop_action_sign_beta_list)
        << "}";
    if (i + 1 < patterns.size()) {
      oss << ",";
    }
  }
  oss << "]";
  return oss.str();
}

}  // namespace

namespace detail {

unsigned resolve_worker_count(const SweepParameterConfig& config) {
  unsigned hw = std::max(1u, std::thread::hardware_concurrency());
  unsigned max_workers = std::max(1u, (hw * 3) / 4);
  if (config.max_workers_override > 0) {
    max_workers = config.max_workers_override;
  }
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

  const bool quiet_logs = config.quiet_logs;
  const unsigned max_workers =
      max_workers_hint > 0 ? max_workers_hint : resolve_worker_count(config);
  Semaphore sem(max_workers);

  newcode::Params base_params = config.base_params;
  base_params.BITGEN_RANDOM_BITS = config.generate_random_bits;
  base_params.NORMALIZE_KNOWN_PREFIX_TAIL = config.normalize_known_prefix_tail;
  newcode::PipelineConfig pipeline_cfg = make_pipeline_config(config);

  std::vector<std::size_t> runnable_indices;
  runnable_indices.reserve(scenarios.size());
  for (std::size_t idx = 0; idx < scenarios.size(); ++idx) {
    const auto& scenario = scenarios[idx];
    if (scenario.alpha_list.size() != base_params.TILES_PER_WIN ||
        scenario.beta_list.size() != base_params.TILES_PER_WIN ||
        scenario.early_stop_action_sign_beta_list.size() != base_params.TILES_PER_WIN) {
      if (log) {
        *log << "[WARN] Scenario '" << scenario.name
             << "' skipped due to list size mismatch (expected "
             << base_params.TILES_PER_WIN << ")\n";
      }
      continue;
    }
    runnable_indices.push_back(idx);
  }

  if (runnable_indices.empty()) {
    if (log) {
      *log << "[ERROR] No scenarios submitted for execution.\n";
    }
    return outputs;
  }

  const auto start_time = std::chrono::steady_clock::now();
  const std::size_t total_jobs = runnable_indices.size();
  std::atomic<std::size_t> finished{0};
  std::mutex progress_mutex;

  const std::string progress_prefix =
      stage_tag.empty() ? "[PROGRESS]" : ("[PROGRESS][" + stage_tag + "]");

  std::vector<std::future<ScenarioOutput>> futures;
  futures.reserve(runnable_indices.size());

  for (std::size_t idx : runnable_indices) {
    const auto& scenario = scenarios[idx];
    newcode::Params params = base_params;
    params.ALPHA_LIST = scenario.alpha_list;
    params.beta_list = scenario.beta_list;
    params.EARLY_STOP_ACTION_SIGN_BETA_LIST =
        scenario.early_stop_action_sign_beta_list;
    if (!params.ALPHA_LIST.empty()) {
      params.ALPHA = params.ALPHA_LIST.front();
    }
    if (!params.beta_list.empty()) {
      params.beta = params.beta_list.front();
    }
    if (!params.EARLY_STOP_ACTION_SIGN_BETA_LIST.empty()) {
      params.EARLY_STOP_ACTION_SIGN_BETA =
          params.EARLY_STOP_ACTION_SIGN_BETA_LIST.front();
    }
    params.CHASE_L = scenario.chase_L;
    params.CHASE_NTEST = scenario.chase_n_test;
    params.CHASE_TOPK_KEEP = scenario.chase_topk_keep;
    params.CHASE_GROUP_MINIMA_BITS = scenario.chase_group_minima_bits;
    params.EARLY_STOP_CONDITION_MODE = scenario.early_stop_condition_mode;
    params.EARLY_STOP_ACTION_MODE = scenario.early_stop_action_mode;
    params.EARLY_STOP_COND_V1_REQUIRE_BCH =
        scenario.early_stop_cond_v1_require_bch;
    params.EARLY_STOP_COND_V1_REQUIRE_OVERALL =
        scenario.early_stop_cond_v1_require_overall;
    params.EARLY_STOP_V2_LLR_ABS_THRESHOLD =
        scenario.early_stop_v2_llr_abs_threshold;
    params.EARLY_STOP_V2_MAX_UNRELIABLE_BITS =
        scenario.early_stop_v2_max_unreliable_bits;
    params.EARLY_STOP_COND_V2_INCLUDE_OVERALL =
        scenario.early_stop_cond_v2_include_overall;
    params.EARLY_STOP_ACTION_HARD_LLR_MAG =
        scenario.early_stop_action_hard_llr_mag;
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
         stage_tag,
         log,
         start_time,
         &finished,
         total_jobs,
         &progress_mutex,
         quiet_logs,
         progress_prefix]() -> ScenarioOutput {
          struct Releaser {
            detail::Semaphore& sem_ref;
            ~Releaser() { sem_ref.release(); }
          } releaser{sem};

          ScenarioOutput output;
          output.idx = idx;
          output.name = scenario.name;
          output.decoder_name = scenario.decoder_name;
          output.alpha_list = scenario.alpha_list;
          output.beta_list = scenario.beta_list;
          output.alpha_start = scenario.alpha_start;
          output.alpha_step = scenario.alpha_step;
          output.beta_start = scenario.beta_start;
          output.beta_step = scenario.beta_step;
          output.early_stop_action_sign_beta_list =
              scenario.early_stop_action_sign_beta_list;
          output.early_stop_beta_start = scenario.early_stop_beta_start;
          output.early_stop_beta_step = scenario.early_stop_beta_step;
          output.early_stop_action_hard_llr_mag =
              scenario.early_stop_action_hard_llr_mag;
          output.chase_L = scenario.chase_L;
          output.chase_n_test = scenario.chase_n_test;
          output.chase_topk_keep = scenario.chase_topk_keep;
          output.chase_group_minima_bits = scenario.chase_group_minima_bits;
          output.mux_group_g = scenario.mux_group_g;
          output.mux_bypass_scheme = scenario.mux_bypass_scheme;
          output.early_stop_condition_mode = scenario.early_stop_condition_mode;
          output.early_stop_action_mode = scenario.early_stop_action_mode;
          output.early_stop_cond_v1_require_bch =
              scenario.early_stop_cond_v1_require_bch;
          output.early_stop_cond_v1_require_overall =
              scenario.early_stop_cond_v1_require_overall;
          output.early_stop_v2_llr_abs_threshold =
              scenario.early_stop_v2_llr_abs_threshold;
          output.early_stop_v2_max_unreliable_bits =
              scenario.early_stop_v2_max_unreliable_bits;
          output.early_stop_cond_v2_include_overall =
              scenario.early_stop_cond_v2_include_overall;
          output.bitgen_seed = scenario.bitgen_seed;
          output.channel_seed = scenario.channel_seed;
          output.ebn0_db = scenario.ebn0_db;

          newcode::Params local_params = params;
          newcode::PipelineConfig local_cfg = pipeline_cfg;
          local_cfg.decoder_name = scenario.decoder_name;
          const std::string label =
              stage_tag.empty() ? scenario.name : (scenario.name + "_" + stage_tag);
          output.result =
              newcode::run_pipeline(local_params, local_cfg, label, scenario.ebn0_db);
          output.ebn0_db = output.result.ebn0_db;
          const auto done = finished.fetch_add(1) + 1;
          const bool need_console_line = quiet_logs;
          const bool need_log_line = (log != nullptr);
          if (need_console_line || need_log_line) {
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
            }
            if (quiet_logs) {
              std::cout << line << "\n";
            }
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
      }
    } catch (...) {
      if (log) {
        *log << "[ERROR] worker threw unknown exception\n";
      }
    }
  }

  if (outputs.empty()) {
    if (log) {
      *log << "[ERROR] No scenario completed successfully.\n";
    }
    return outputs;
  }

  std::sort(outputs.begin(), outputs.end(),
            [](const auto& a, const auto& b) { return a.idx < b.idx; });

  return outputs;
}

}  // namespace detail

int run_sweep(const SweepParameterConfig& config) {
  using namespace detail;
  SweepParameterConfig resolved = config;
  if (!resolved.siso_active_list.empty()) {
    resolved.base_params.SISO_ACTIVE_LIST = resolved.siso_active_list;
  }
  resolved.base_params.ENABLE_EARLY_STOP = resolved.enable_early_stop;
  resolved.base_params.EARLY_STOP_CONDITION_MODE =
      resolved.early_stop_condition_mode;
  resolved.base_params.EARLY_STOP_ACTION_MODE =
      resolved.early_stop_action_mode;
  resolved.base_params.EARLY_STOP_COND_V1_REQUIRE_BCH =
      resolved.early_stop_cond_v1_require_bch;
  resolved.base_params.EARLY_STOP_COND_V1_REQUIRE_OVERALL =
      resolved.early_stop_cond_v1_require_overall;
  resolved.base_params.EARLY_STOP_V2_LLR_ABS_THRESHOLD =
      resolved.early_stop_v2_llr_abs_threshold;
  resolved.base_params.EARLY_STOP_V2_MAX_UNRELIABLE_BITS =
      resolved.early_stop_v2_max_unreliable_bits;
  resolved.base_params.EARLY_STOP_COND_V2_INCLUDE_OVERALL =
      resolved.early_stop_cond_v2_include_overall;
  resolved.base_params.EARLY_STOP_ACTION_RESIDUAL_DIVISOR =
      resolved.early_stop_action_residual_divisor;
  resolved.base_params.EARLY_STOP_ACTION_HARD_LLR_MAG =
      resolved.early_stop_action_hard_llr_mag;
  resolved.base_params.CHASE_NTEST = resolved.chase_n_test;
  resolved.base_params.CHASE_TOPK_KEEP = resolved.chase_topk_keep;
  resolved.base_params.CHASE_GROUP_MINIMA_BITS =
      resolved.chase_group_minima_bits;
  resolved.base_params.MUX_GROUP_G = resolved.mux_group_g;
  resolved.base_params.MUX_ENABLE_RECONFIG = resolved.mux_enable_reconfig;
  resolved.base_params.MUX_EXTRA_BYPASS_EDGES =
      resolved.mux_extra_bypass_edges;
  if (std::isfinite(resolved.early_stop_action_sign_beta_fill)) {
    resolved.base_params.EARLY_STOP_ACTION_SIGN_BETA =
        resolved.early_stop_action_sign_beta_fill;
  }
  auto validate_condition_mode = [](int mode) -> bool {
    return mode == 1 || mode == 2;
  };
  auto validate_action_mode = [](int mode) -> bool {
    return mode == 1 || mode == 2 || mode == 3;
  };
  const int resolved_condition_mode =
      resolved.base_params.EARLY_STOP_CONDITION_MODE;
  if (!validate_condition_mode(resolved_condition_mode)) {
    std::cerr << "[ERROR] early_stop_condition_mode must be 1 or 2\n";
    return 1;
  }
  if (!validate_action_mode(resolved.early_stop_action_mode)) {
    std::cerr << "[ERROR] early_stop_action_mode must be 1, 2 or 3\n";
    return 1;
  }
  for (int mode : resolved.early_stop_condition_candidates) {
    if (!validate_condition_mode(mode)) {
      std::cerr << "[ERROR] early_stop_condition_candidates contains invalid mode\n";
      return 1;
    }
  }
  for (int mode : resolved.early_stop_action_candidates) {
    if (!validate_action_mode(mode)) {
      std::cerr << "[ERROR] early_stop_action_candidates contains invalid mode\n";
      return 1;
    }
  }
  if (resolved.early_stop_v2_llr_abs_threshold < 0.0f) {
    std::cerr << "[ERROR] early_stop_v2_llr_abs_threshold must be >= 0\n";
    return 1;
  }
  for (float threshold : resolved.early_stop_v2_llr_abs_threshold_candidates) {
    if (threshold < 0.0f) {
      std::cerr << "[ERROR] early_stop_v2_llr_abs_threshold_candidates must be >= 0\n";
      return 1;
    }
  }
  if (resolved.early_stop_v2_max_unreliable_bits < 0 ||
      resolved.early_stop_v2_max_unreliable_bits >
          static_cast<int>(newcode::Params::BCH_N)) {
    std::cerr << "[ERROR] early_stop_v2_max_unreliable_bits must be in [0, BCH_N]\n";
    return 1;
  }
  for (int count : resolved.early_stop_v2_max_unreliable_bits_candidates) {
    if (count < 0 ||
        count > static_cast<int>(newcode::Params::BCH_N)) {
      std::cerr << "[ERROR] early_stop_v2_max_unreliable_bits_candidates must be in [0, BCH_N]\n";
      return 1;
    }
  }
  if (!(resolved.early_stop_action_residual_divisor > 0.0f)) {
    std::cerr << "[ERROR] early_stop_action_residual_divisor must be > 0\n";
    return 1;
  }
  if (!std::isfinite(resolved.early_stop_action_hard_llr_mag)) {
    std::cerr << "[ERROR] early_stop_action_hard_llr_mag must be finite\n";
    return 1;
  }
  if (resolved.chase_topk_keep < 1) {
    std::cerr << "[ERROR] chase_topk_keep must be >= 1\n";
    return 1;
  }
  if (resolved.chase_n_test < 1) {
    std::cerr << "[ERROR] chase_n_test must be >= 1\n";
    return 1;
  }
  if (resolved.chase_group_minima_bits < 0) {
    std::cerr << "[ERROR] chase_group_minima_bits must be >= 0\n";
    return 1;
  }
  for (float hard_mag : resolved.early_stop_action_hard_llr_mag_candidates) {
    if (!std::isfinite(hard_mag)) {
      std::cerr << "[ERROR] early_stop_action_hard_llr_mag_candidates must be finite\n";
      return 1;
    }
  }
  for (int keep : resolved.chase_topk_keep_candidates) {
    if (keep < 1) {
      std::cerr << "[ERROR] chase_topk_keep_candidates must be >= 1\n";
      return 1;
    }
  }
  for (int n_test : resolved.chase_n_test_candidates) {
    if (n_test < 1) {
      std::cerr << "[ERROR] chase_n_test_candidates must be >= 1\n";
      return 1;
    }
  }
  for (int group_bits : resolved.chase_group_minima_bits_candidates) {
    if (group_bits < 0) {
      std::cerr << "[ERROR] chase_group_minima_bits_candidates must be >= 0\n";
      return 1;
    }
  }
  const auto mux_ok = newcode::mux::validate_siso_active_list(
      resolved.base_params.SISO_ACTIVE_LIST,
      resolved.base_params.TILES_PER_WIN);
  if (!mux_ok.ok) {
    std::cerr << "[ERROR] " << mux_ok.error << "\n";
    return 1;
  }
  const std::size_t rows_to_decode =
      static_cast<std::size_t>(resolved.base_params.CHASE_SBR) *
      newcode::Params::BITS_PER_SUBBLOCK_DIM;
  const auto group_ok = newcode::mux::validate_mux_group_g(
      resolved.base_params.MUX_GROUP_G,
      rows_to_decode);
  if (!group_ok.ok) {
    std::cerr << "[ERROR] " << group_ok.error << "\n";
    return 1;
  }
  if (resolved.base_params.MUX_ENABLE_RECONFIG) {
    for (std::size_t tile_idx = 0;
         tile_idx < resolved.base_params.SISO_ACTIVE_LIST.size();
         ++tile_idx) {
      const auto reconfig_ok = newcode::mux::validate_mux_reconfig_runtime(
          resolved.base_params.MUX_GROUP_G,
          resolved.base_params.SISO_ACTIVE_LIST[tile_idx],
          rows_to_decode);
      if (!reconfig_ok.ok) {
        std::cerr << "[ERROR] tile " << tile_idx << ": "
                  << reconfig_ok.error << "\n";
        return 1;
      }
    }
  }
  const SweepParameterConfig& cfg = resolved;

  const std::filesystem::path data_dir = "data";
  io::ensure_dir(data_dir);
  const std::string run_id = utils::now_stamp();
  const std::string log_path = (data_dir / ("run_" + run_id + ".log")).string();
  const bool mirror_console = !cfg.quiet_logs;
  DualOut out(std::cout, log_path, mirror_console);

  out << "[CONFIG] ofec_sweep top-level parameters\n";
  out << "[CONFIG][TX] interleaver_name=" << cfg.interleaver_name
      << " bits_per_symbol=" << cfg.bits_per_symbol
      << " generate_random_bits=" << (cfg.generate_random_bits ? "true" : "false")
      << " bitgen_seed=" << cfg.base_params.BITGEN_SEED
      << " bitgen_seed_count=" << cfg.bitgen_seed_count
      << " bitgen_seed_candidates=" << join_compact(cfg.bitgen_seed_candidates)
      << "\n";
  out << "[CONFIG][CH] ebn0_start=" << cfg.ebn0_start
      << " ebn0_end=" << cfg.ebn0_end
      << " ebn0_points=" << cfg.ebn0_points
      << " ebn0_candidates=" << join_compact(cfg.ebn0_candidates)
      << " channel_seed=" << cfg.base_params.CHANNEL_SEED
      << " channel_seed_count=" << cfg.channel_seed_count
      << " channel_seed_candidates=" << join_compact(cfg.channel_seed_candidates)
      << "\n";
  out << "[CONFIG][Q] llr_bits=" << cfg.base_params.LLR_BITS
      << " quant_clip_ratio=" << cfg.quant_clip_ratio << "\n";
  out << "[CONFIG][DEC] decoder_name=" << cfg.decoder_name
      << " decoder_name_candidates=" << join_compact(cfg.decoder_name_candidates)
      << " normalize_extrinsic=" << (cfg.normalize_extrinsic ? "true" : "false")
      << " normalize_known_prefix_tail="
      << (cfg.normalize_known_prefix_tail ? "true" : "false")
      << " chase_l_candidates=" << join_compact(cfg.chase_l_candidates)
      << " chase_n_test=" << cfg.chase_n_test
      << " chase_n_test_candidates=" << join_compact(cfg.chase_n_test_candidates)
      << " chase_topk_keep=" << cfg.chase_topk_keep
      << " chase_topk_keep_candidates=" << join_compact(cfg.chase_topk_keep_candidates)
      << " chase_group_minima_bits=" << cfg.chase_group_minima_bits
      << " chase_group_minima_bits_candidates="
      << join_compact(cfg.chase_group_minima_bits_candidates)
      << " alpha_start_candidates=" << join_compact(cfg.alpha_start_candidates)
      << " alpha_step_candidates=" << join_compact(cfg.alpha_step_candidates)
      << " beta_start_candidates=" << join_compact(cfg.beta_start_candidates)
      << " beta_step_candidates=" << join_compact(cfg.beta_step_candidates)
      << " explicit_patterns=" << summarize_patterns(cfg.explicit_patterns)
      << " siso_active_list=" << join_compact(cfg.siso_active_list)
      << " mux_group_g=" << cfg.mux_group_g
      << " mux_enable_reconfig=" << (cfg.mux_enable_reconfig ? "true" : "false")
      << " mux_bypass_scheme=" << cfg.mux_bypass_scheme
      << "\n";
  out << "[CONFIG][ES] enable_early_stop=" << (cfg.enable_early_stop ? "true" : "false")
      << " early_stop_condition_mode=" << cfg.early_stop_condition_mode
      << " early_stop_condition_candidates="
      << join_compact(cfg.early_stop_condition_candidates)
      << " early_stop_action_mode=" << cfg.early_stop_action_mode
      << " early_stop_action_candidates=" << join_compact(cfg.early_stop_action_candidates)
      << " v1_require_bch=" << (cfg.early_stop_cond_v1_require_bch ? "true" : "false")
      << " v1_require_bch_candidates="
      << join_compact(cfg.early_stop_cond_v1_require_bch_candidates)
      << " v1_require_overall=" << (cfg.early_stop_cond_v1_require_overall ? "true" : "false")
      << " v1_require_overall_candidates="
      << join_compact(cfg.early_stop_cond_v1_require_overall_candidates)
      << " v2_llr_abs_threshold=" << cfg.early_stop_v2_llr_abs_threshold
      << " v2_llr_abs_threshold_candidates="
      << join_compact(cfg.early_stop_v2_llr_abs_threshold_candidates)
      << " v2_max_unreliable_bits=" << cfg.early_stop_v2_max_unreliable_bits
      << " v2_max_unreliable_bits_candidates="
      << join_compact(cfg.early_stop_v2_max_unreliable_bits_candidates)
      << " v2_include_overall=" << (cfg.early_stop_cond_v2_include_overall ? "true" : "false")
      << " action_sign_beta_fill=" << cfg.early_stop_action_sign_beta_fill
      << " action_beta_start_candidates="
      << join_compact(cfg.early_stop_action_beta_start_candidates)
      << " action_beta_step_candidates="
      << join_compact(cfg.early_stop_action_beta_step_candidates)
      << " action_residual_divisor=" << cfg.early_stop_action_residual_divisor
      << " action_hard_llr_mag=" << cfg.early_stop_action_hard_llr_mag
      << " action_hard_llr_mag_candidates="
      << join_compact(cfg.early_stop_action_hard_llr_mag_candidates)
      << "\n";
  out << "[CONFIG][DBG] quiet_pipeline=" << (cfg.quiet_pipeline ? "true" : "false")
      << " quiet_logs=" << (cfg.quiet_logs ? "true" : "false")
      << " trace_enable=" << (cfg.base_params.debug_trace.enable ? "true" : "false")
      << " trace_row=" << cfg.base_params.debug_trace.row
      << " trace_col=" << cfg.base_params.debug_trace.col
      << " trace_log_read=" << (cfg.base_params.debug_trace.log_read_mapping ? "true" : "false")
      << " trace_log_write=" << (cfg.base_params.debug_trace.log_write_mapping ? "true" : "false")
      << " trace_log_mismatch=" << (cfg.base_params.debug_trace.log_mismatch ? "true" : "false")
      << "\n";
  out << std::flush;

  const std::string csv_path =
      (data_dir / ("ofec_sweep_results_" + run_id + ".csv")).string();
  ensure_csv_header(csv_path);

  std::vector<int> bitgen_seeds =
      resolve_seeds(cfg.bitgen_seed_candidates,
                    cfg.bitgen_seed_count,
                    cfg.base_params.BITGEN_SEED);
  std::vector<int> channel_seeds =
      resolve_seeds(cfg.channel_seed_candidates,
                    cfg.channel_seed_count,
                    cfg.base_params.CHANNEL_SEED);
  const std::vector<float> ebn0_values = build_ebn0_values(cfg);

  auto scenarios = build_scenarios(cfg, ebn0_values, bitgen_seeds, channel_seeds);
  const std::size_t scenario_count = scenarios.size();
  std::cout << "[INFO] total scenarios = " << scenario_count << "\n";
  if (scenario_count == 0) {
    out << "[ERROR] No scenarios generated.\n";
    return 1;
  }

  const unsigned max_workers = resolve_worker_count(cfg);
  std::cout << "[INFO] using up to " << max_workers << " workers\n";
  auto results = run_scenarios_parallel(scenarios, cfg, max_workers, "", &out);
  if (results.empty()) {
    return 1;
  }

  std::ofstream csv(csv_path, std::ios::out | std::ios::app);
  csv.setf(std::ios::fixed);
  csv << std::setprecision(8);

  std::vector<std::string> scenario_summaries;
  scenario_summaries.reserve(results.size());

  double best_post_ber = std::numeric_limits<double>::infinity();
  std::size_t best_index = static_cast<std::size_t>(-1);
  newcode::PipelineResult best_result{};
  std::vector<double> best_tile_early_stop_pct;
  std::vector<double> best_tile_row_early_stop_pct;
  float best_alpha_start = 0.0f;
  float best_alpha_step = 0.0f;
  float best_beta_start = 0.0f;
  float best_beta_step = 0.0f;
  int best_chase_L = cfg.base_params.CHASE_L;
  int best_chase_n_test = cfg.base_params.CHASE_NTEST;
  int best_chase_topk_keep = cfg.base_params.CHASE_TOPK_KEEP;
  int best_chase_group_minima_bits = cfg.base_params.CHASE_GROUP_MINIMA_BITS;
  int best_bitgen_seed = cfg.base_params.BITGEN_SEED;
  int best_channel_seed = cfg.base_params.CHANNEL_SEED;
  float best_ebn0_db = !ebn0_values.empty() ? ebn0_values.front() : newcode::DEFAULT_EBN0_DB;

  for (const auto& pack : results) {
    const auto& result = pack.result;
    if (!result.tile_early_stop_pct.empty()) {
      out << "[INFO] " << pack.name << " tile early-stop hit rates (%): ";
      out << std::fixed << std::setprecision(1);
      for (size_t i = 0; i < result.tile_early_stop_pct.size(); ++i) {
        out << result.tile_early_stop_pct[i]
            << (i + 1 < result.tile_early_stop_pct.size() ? ", " : "\n");
      }
      out << std::defaultfloat;
    }
    if (!result.tile_row_early_stop_pct.empty()) {
      out << "[INFO] " << pack.name << " tile row-early-stop hit rates (%): ";
      out << std::fixed << std::setprecision(1);
      for (size_t i = 0; i < result.tile_row_early_stop_pct.size(); ++i) {
        out << result.tile_row_early_stop_pct[i]
            << (i + 1 < result.tile_row_early_stop_pct.size() ? ", " : "\n");
      }
      out << std::defaultfloat;
    }

    std::ostringstream summary;
    summary << "[SUMMARY] " << pack.name
            << " decoder=" << pack.decoder_name
            << " Pre-FEC BER=" << result.pre_fec.ber
            << " (errs=" << result.pre_fec.errors << "/" << result.pre_fec.total << ")"
            << " | Post-FEC BER=" << result.post_fec.ber
            << " (errs=" << result.post_fec.errors << "/" << result.post_fec.total << ")";
    summary << std::fixed << std::setprecision(3)
            << " | alpha_start=" << pack.alpha_start
            << " alpha_step=" << pack.alpha_step
            << " | beta_start=" << pack.beta_start
            << " beta_step=" << pack.beta_step
            << " | es_beta_start=" << pack.early_stop_beta_start
            << " es_beta_step=" << pack.early_stop_beta_step
            << " | hard_mag=" << pack.early_stop_action_hard_llr_mag
            << " | Eb/N0=" << pack.ebn0_db
            << std::defaultfloat
            << " | CHASE_L=" << pack.chase_L
            << " CHASE_NTEST=" << pack.chase_n_test
            << " CHASE_TOPK_KEEP=" << pack.chase_topk_keep
            << " CHASE_GROUP_MINIMA_BITS=" << pack.chase_group_minima_bits
            << " | cond/action=" << pack.early_stop_condition_mode
            << "/" << pack.early_stop_action_mode
            << " | Seeds(bit/channel)=" << pack.bitgen_seed << "/" << pack.channel_seed;
    if (!result.tile_early_stop_pct.empty()) {
      summary << " | EarlyStop%=["
              << detail::join_vec(result.tile_early_stop_pct, ',', 1)
              << "]";
    }
    if (!result.tile_row_early_stop_pct.empty()) {
      summary << " | EarlyStopRow%=["
              << detail::join_vec(result.tile_row_early_stop_pct, ',', 1)
              << "]";
    }
    scenario_summaries.push_back(summary.str());
    out << summary.str() << "\n";

    write_csv_row(csv,
                  utils::now_stamp(),
                  run_id,
                  "" /*stage*/,
                  0 /*num_bits*/,
                  scenarios[pack.idx],
                  result,
                  CsvFormat::Basic);

    if (result.post_fec.total > 0 && result.post_fec.ber < best_post_ber) {
      best_post_ber = result.post_fec.ber;
      best_index = pack.idx;
      best_result = result;
      best_tile_early_stop_pct = result.tile_early_stop_pct;
      best_tile_row_early_stop_pct = result.tile_row_early_stop_pct;
      best_alpha_start = pack.alpha_start;
      best_alpha_step = pack.alpha_step;
      best_beta_start = pack.beta_start;
      best_beta_step = pack.beta_step;
      best_chase_L = pack.chase_L;
      best_chase_n_test = pack.chase_n_test;
      best_chase_topk_keep = pack.chase_topk_keep;
      best_chase_group_minima_bits = pack.chase_group_minima_bits;
      best_bitgen_seed = pack.bitgen_seed;
      best_channel_seed = pack.channel_seed;
      best_ebn0_db = pack.ebn0_db;
    }
  }

  if (best_index == static_cast<std::size_t>(-1)) {
    out << "[ERROR] No valid scenarios evaluated.\n";
    return 1;
  }

  out << "\n[SUMMARY] All scenarios:\n";
  for (const auto& line : scenario_summaries) {
    out << "  " << line << '\n';
  }

  const auto& best_scenario = scenarios[best_index];
  out << "\n[RESULT] Best scenario: " << best_scenario.name
      << " decoder=" << best_scenario.decoder_name
      << " with Post-FEC BER=" << best_result.post_fec.ber
      << " (errs=" << best_result.post_fec.errors << "/" << best_result.post_fec.total << ")\n";
  out << std::fixed << std::setprecision(3);
  out << "[RESULT] Best alpha start/step: start=" << best_alpha_start
      << " step=" << best_alpha_step << "\n";
  out << "[RESULT] Best beta start/step: start=" << best_beta_start
      << " step=" << best_beta_step << "\n";
  out << "[RESULT] Best sweep Eb/N0: " << best_ebn0_db << " dB\n";
  out << std::defaultfloat;
  out << "[RESULT] Best early-stop cond/action: "
      << best_scenario.early_stop_condition_mode
      << " / " << best_scenario.early_stop_action_mode << "\n";
  if (best_scenario.early_stop_action_mode == 3) {
    out << std::fixed << std::setprecision(3)
        << "[RESULT] Best early-stop hard_mag: "
        << best_scenario.early_stop_action_hard_llr_mag << "\n"
        << std::defaultfloat;
  }
  out << "[RESULT] Best CHASE_L/CHASE_NTEST: " << best_chase_L
      << " / " << best_chase_n_test << "\n";
  out << "[RESULT] Best CHASE_TOPK_KEEP: " << best_chase_topk_keep << "\n";
  out << "[RESULT] Best CHASE_GROUP_MINIMA_BITS: "
      << best_chase_group_minima_bits << "\n";
  out << "[RESULT] Best RNG seeds (bit/channel): "
      << best_bitgen_seed << "/" << best_channel_seed << "\n";

  if (!best_tile_early_stop_pct.empty()) {
    out << "[RESULT] Best tile early-stop hit rates (%): ";
    out << std::fixed << std::setprecision(1);
    for (size_t i = 0; i < best_tile_early_stop_pct.size(); ++i) {
      out << best_tile_early_stop_pct[i]
          << (i + 1 < best_tile_early_stop_pct.size() ? ", " : "\n");
    }
    out << std::defaultfloat;
  }
  if (!best_tile_row_early_stop_pct.empty()) {
    out << "[RESULT] Best tile row-early-stop hit rates (%): ";
    out << std::fixed << std::setprecision(1);
    for (size_t i = 0; i < best_tile_row_early_stop_pct.size(); ++i) {
      out << best_tile_row_early_stop_pct[i]
          << (i + 1 < best_tile_row_early_stop_pct.size() ? ", " : "\n");
    }
    out << std::defaultfloat;
  }

  out << "[RESULT] Best ALPHA_LIST: ";
  for (size_t i = 0; i < best_scenario.alpha_list.size(); ++i) {
    out << best_scenario.alpha_list[i]
        << (i + 1 < best_scenario.alpha_list.size() ? ", " : "\n");
  }

  out << "[RESULT] Best beta_list: ";
  for (size_t i = 0; i < best_scenario.beta_list.size(); ++i) {
    out << best_scenario.beta_list[i]
        << (i + 1 < best_scenario.beta_list.size() ? ", " : "\n");
  }

  out << "[RESULT] Best early-stop beta_list: ";
  for (size_t i = 0; i < best_scenario.early_stop_action_sign_beta_list.size(); ++i) {
    out << best_scenario.early_stop_action_sign_beta_list[i]
        << (i + 1 < best_scenario.early_stop_action_sign_beta_list.size() ? ", " : "\n");
  }

  return 0;
}

}  // namespace ofec_sweep
