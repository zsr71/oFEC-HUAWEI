#pragma once

#include <condition_variable>
#include <fstream>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

#include "newcode/ofec_sweep_runner.hpp"

namespace ofec_sweep {
namespace detail {

struct SweepScenario {
  std::string name;
  std::string decoder_name;
  // ofec_sweep3 使用的 TwoMain 策略维度。其他 sweep 入口不展开该维度，
  // 因而保持 Params 的默认 Legacy 语义。
  std::string twomain_policy_label = "legacy";
  newcode::TwoMainHisoOutputMode twomain_hiso_output_mode =
      newcode::TwoMainHisoOutputMode::Legacy;
  float twomain_hiso_m2 = 99.0f;
  float twomain_hiso_rho_corr = 1.0f;
  float twomain_hiso_rho_keep = 1.0f;
  unsigned level56_hiso_allowed_class_mask = 0x0fu;
  int level56_shared_hiso_active = 32;
  int level56_shared_siso_active = 32;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  std::vector<float> early_stop_action_sign_beta_list;
  float alpha_start = 0.0f;
  float alpha_step = 0.0f;
  float beta_start = 0.0f;
  float beta_step = 0.0f;
  float early_stop_beta_start = 0.0f;
  float early_stop_beta_step = 0.0f;
  float early_stop_action_hard_llr_mag = 1.0f;
  float alpha_low = 0.0f;
  float alpha_high = 0.0f;
  float gamma_alpha = 1.0f;
  float beta_low = 0.0f;
  float beta_high = 0.0f;
  float gamma_beta = 1.0f;
  int chase_L = 0;
  int chase_n_test = 0;
  int chase_topk_keep = 8;
  int chase_group_minima_bits = 3;
  int mux_group_g = 1;
  int mux_scheduling_mode = 0;
  int mux_early_stop_priority_rule = 0;
  int mux_bypass_scheme = 0;
  int early_stop_condition_mode = 1;
  int early_stop_action_mode = 1;
  int early_stop_bind_group_size = 1;
  bool early_stop_cond_v1_require_bch = true;
  bool early_stop_cond_v1_require_overall = true;
  float early_stop_v2_llr_abs_threshold = 0.5f;
  int early_stop_v2_max_unreliable_bits = 8;
  bool early_stop_cond_v2_include_overall = true;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = newcode::DEFAULT_EBN0_DB;
};

struct ScenarioOutput {
  std::size_t idx{};
  std::string name;
  std::string decoder_name;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  std::vector<float> early_stop_action_sign_beta_list;
  float alpha_start = 0.0f;
  float alpha_step = 0.0f;
  float beta_start = 0.0f;
  float beta_step = 0.0f;
  float early_stop_beta_start = 0.0f;
  float early_stop_beta_step = 0.0f;
  float early_stop_action_hard_llr_mag = 1.0f;
  int chase_L = 0;
  int chase_n_test = 0;
  int chase_topk_keep = 8;
  int chase_group_minima_bits = 3;
  int mux_group_g = 1;
  int mux_scheduling_mode = 0;
  int mux_early_stop_priority_rule = 0;
  int mux_bypass_scheme = 0;
  int early_stop_condition_mode = 1;
  int early_stop_action_mode = 1;
  int early_stop_bind_group_size = 1;
  bool early_stop_cond_v1_require_bch = true;
  bool early_stop_cond_v1_require_overall = true;
  float early_stop_v2_llr_abs_threshold = 0.5f;
  int early_stop_v2_max_unreliable_bits = 8;
  bool early_stop_cond_v2_include_overall = true;
  int bitgen_seed = 0;
  int channel_seed = 0;
  float ebn0_db = newcode::DEFAULT_EBN0_DB;
  newcode::PipelineResult result;
};

class DualOut {
 public:
  DualOut(std::ostream& console, const std::string& filepath, bool mirror_console = true);

  template <typename T>
  DualOut& operator<<(const T& value) {
    if (console_) {
      *console_ << value;
    }
    if (file_) {
      file_ << value;
    }
    return *this;
  }

  DualOut& operator<<(std::ostream& (*pf)(std::ostream&));

 private:
  std::ostream* console_;
  std::ofstream file_;
};

class Semaphore {
 public:
  explicit Semaphore(std::size_t count);

  void acquire();
  void release();

 private:
  std::mutex mutex_;
  std::condition_variable cv_;
  std::size_t count_;
};

void ensure_csv_header(const std::string& csv_path);
void ensure_csv_header_v2(const std::string& csv_path);

std::vector<float> generate_sequence(float start, float step, std::size_t length);
float infer_step(const std::vector<float>& values);
std::vector<int> generate_random_seeds(int count);

newcode::PipelineConfig make_pipeline_config(const SweepParameterConfig& config);

std::vector<SweepScenario> build_scenarios(const SweepParameterConfig& config,
                                           const std::vector<float>& ebn0_candidates,
                                           const std::vector<int>& bitgen_seeds,
                                           const std::vector<int>& channel_seeds);

std::vector<float> build_ebn0_values(const SweepParameterConfig& config);

std::string join_vec(const std::vector<float>& values, char sep, int precision);
std::string join_vec(const std::vector<double>& values, char sep, int precision);
double mean(const std::vector<double>& values);
enum class CsvFormat {
  Basic,
  Extended
};

struct ExtendedCsvConfigSnapshot {
  std::string bitgen_seed_candidates;
  std::string channel_seed_candidates;
  std::string ebn0_candidates;
  std::string decoder_name_candidates;
  std::string chase_l_candidates;
  std::string chase_n_test_candidates;
  std::string chase_topk_keep_candidates;
  std::string chase_group_minima_bits_candidates;
  std::string alpha_start_candidates;
  std::string alpha_step_candidates;
  std::string beta_start_candidates;
  std::string beta_step_candidates;
  std::string explicit_patterns;
  std::string early_stop_condition_candidates;
  std::string early_stop_condition_mode_list;
  std::string early_stop_action_candidates;
  std::string early_stop_action_mode_list;
  std::string early_stop_bind_group_size_list;
  std::string early_stop_cond_v1_require_bch_candidates;
  std::string early_stop_cond_v1_require_overall_candidates;
  std::string early_stop_v2_llr_abs_threshold_candidates;
  std::string early_stop_v2_max_unreliable_bits_candidates;
  std::string early_stop_action_beta_start_candidates;
  std::string early_stop_action_beta_step_candidates;
  std::string early_stop_action_hard_llr_mag_candidates;
  std::size_t tiles_per_window = 0;
  float eval_ebn0_db = 0.0f;
  std::size_t stage1_bits = 0;
  std::size_t stage2_bits = 0;
  float keep_ratio = 0.0f;
  std::string interleaver_name;
  int bitgen_seed = 0;
  int bitgen_seed_count = 0;
  float ebn0_start = 0.0f;
  float ebn0_end = 0.0f;
  int ebn0_points = 0;
  int channel_seed = 0;
  int channel_seed_count = 0;
  unsigned bits_per_symbol = 0;
  bool normalize_extrinsic = false;
  bool generate_random_bits = false;
  bool normalize_known_prefix_tail = false;
  std::size_t llr_bits = 0;
  float quant_clip_ratio = 0.0f;
  int chase_n_test = 0;
  int chase_topk_keep = 0;
  int chase_group_minima_bits = 0;
  std::vector<int> siso_active_list;
  int mux_group_g = 0;
  int mux_scheduling_mode = 0;
  std::string mux_scheduling_mode_candidates;
  int mux_early_stop_priority_rule = 0;
  std::string mux_early_stop_priority_rule_candidates;
  bool mux_enable_reconfig = false;
  int mux_bypass_scheme = 0;
  bool enable_early_stop = false;
  int early_stop_bind_group_size = 1;
  int early_stop_condition_mode = 0;
  int early_stop_action_mode = 0;
  bool early_stop_cond_v1_require_bch = false;
  bool early_stop_cond_v1_require_overall = false;
  float early_stop_v2_llr_abs_threshold = 0.0f;
  int early_stop_v2_max_unreliable_bits = 0;
  bool early_stop_cond_v2_include_overall = false;
  float early_stop_action_sign_beta_fill = 0.0f;
  float early_stop_action_residual_divisor = 0.0f;
  float early_stop_action_hard_llr_mag = 0.0f;
  bool quiet_pipeline = false;
  bool quiet_logs = false;
  bool trace_enable = false;
  long trace_row = -1;
  long trace_col = -1;
  bool trace_log_read = false;
  bool trace_log_write = false;
  bool trace_log_mismatch = false;
  std::vector<float> alpha_low_grid;
  std::vector<float> alpha_high_grid;
  std::vector<float> beta_low_grid;
  std::vector<float> beta_high_grid;
  std::vector<float> gamma_alpha_grid;
  std::vector<float> gamma_beta_grid;
};

void write_csv_row(std::ostream& csv,
                   const std::string& timestamp,
                   const std::string& run_id,
                   const std::string& stage_tag,
                   std::size_t num_bits,
                   const SweepScenario& scenario,
                   const newcode::PipelineResult& result,
                   CsvFormat format,
                   const ExtendedCsvConfigSnapshot* snapshot = nullptr);

std::vector<ScenarioOutput> run_scenarios_parallel(
    const std::vector<SweepScenario>& scenarios,
    const SweepParameterConfig& config,
    unsigned max_workers_hint = 0,
    const std::string& stage_tag = std::string{},
    DualOut* log = nullptr);

unsigned resolve_worker_count(const SweepParameterConfig& config);

}  // namespace detail
}  // namespace ofec_sweep
