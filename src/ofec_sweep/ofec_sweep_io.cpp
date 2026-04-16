#include "ofec_sweep_detail.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace ofec_sweep {
namespace detail {

namespace {

std::string join_vec_int(const std::vector<int>& values, char sep) {
  std::ostringstream oss;
  for (size_t i = 0; i < values.size(); ++i) {
    oss << values[i];
    if (i + 1 < values.size()) {
      oss << sep;
    }
  }
  return oss.str();
}

std::string join_vec_size_t(const std::vector<std::size_t>& values, char sep) {
  std::ostringstream oss;
  for (size_t i = 0; i < values.size(); ++i) {
    oss << values[i];
    if (i + 1 < values.size()) {
      oss << sep;
    }
  }
  return oss.str();
}

}  // namespace

DualOut::DualOut(std::ostream& console, const std::string& filepath, bool mirror_console)
    : console_(mirror_console ? &console : nullptr),
      file_(filepath, std::ios::out | std::ios::app) {
  if (console_) {
    *console_ << std::unitbuf;
  }
  if (file_) {
    file_ << std::unitbuf;
  }
}

DualOut& DualOut::operator<<(std::ostream& (*pf)(std::ostream&)) {
  if (console_) {
    pf(*console_);
  }
  if (file_) {
    pf(file_);
  }
  return *this;
}

void ensure_csv_header(const std::string& csv_path) {
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,scenario,decoder_name,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "early_stop_action_hard_llr_mag,"
          "chase_L,chase_n_test,chase_topk_keep,chase_group_minima_bits,mux_group_g,mux_scheduling_mode,mux_early_stop_priority_rule,mux_bypass_scheme,bitgen_seed,channel_seed,ebn0_db,ALPHA_LIST,beta_list,"
          "early_stop_beta_list,early_stop_condition_mode,early_stop_action_mode,"
          "early_stop_bind_group_size,"
          "early_stop_cond_v1_require_bch,early_stop_cond_v1_require_overall,"
          "early_stop_v2_llr_abs_threshold,early_stop_v2_max_unreliable_bits,"
          "early_stop_cond_v2_include_overall,"
          "early_stop_mean_pct,early_stop_row_mean_pct,early_stop_list,early_stop_row_list,"
          "unscheduled_mean_pct,unscheduled_need_mean_pct,unscheduled_count_list,unscheduled_list,unscheduled_need_list,"
          "interleaver_name,bits_per_symbol,generate_random_bits,bitgen_seed_cfg,bitgen_seed_count,bitgen_seed_candidates,"
          "ebn0_start,ebn0_end,ebn0_points,ebn0_candidates,channel_seed_cfg,channel_seed_count,channel_seed_candidates,"
          "llr_bits,quant_clip_ratio,normalize_extrinsic,normalize_known_prefix_tail,"
          "decoder_name_candidates,chase_l_candidates,chase_n_test_cfg,chase_n_test_candidates,"
          "chase_topk_keep_cfg,chase_topk_keep_candidates,chase_group_minima_bits_cfg,chase_group_minima_bits_candidates,"
          "alpha_start_candidates,alpha_step_candidates,beta_start_candidates,beta_step_candidates,"
          "explicit_patterns,siso_active_list,mux_group_g_cfg,mux_scheduling_mode_cfg,mux_scheduling_mode_candidates,mux_early_stop_priority_rule_cfg,mux_early_stop_priority_rule_candidates,mux_enable_reconfig,mux_bypass_scheme_cfg,"
          "enable_early_stop,early_stop_condition_mode_cfg,early_stop_condition_candidates,"
          "early_stop_condition_mode_list_cfg,"
          "early_stop_action_mode_cfg,early_stop_action_mode_list_cfg,"
          "early_stop_bind_group_size_cfg,early_stop_bind_group_size_list_cfg,early_stop_action_candidates,"
          "early_stop_cond_v1_require_bch_cfg,early_stop_cond_v1_require_bch_candidates,"
          "early_stop_cond_v1_require_overall_cfg,early_stop_cond_v1_require_overall_candidates,"
          "early_stop_v2_llr_abs_threshold_cfg,early_stop_v2_llr_abs_threshold_candidates,"
          "early_stop_v2_max_unreliable_bits_cfg,early_stop_v2_max_unreliable_bits_candidates,"
          "early_stop_cond_v2_include_overall_cfg,early_stop_action_sign_beta_fill,"
          "early_stop_action_beta_start_candidates,early_stop_action_beta_step_candidates,"
          "early_stop_action_residual_divisor,early_stop_action_hard_llr_mag_cfg,early_stop_action_hard_llr_mag_candidates,"
          "quiet_pipeline,quiet_logs,trace_enable,trace_row,trace_col,trace_log_read,trace_log_write,trace_log_mismatch\n";
}

void ensure_csv_header_v2(const std::string& csv_path) {
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,stage,num_bits,label,"
          "tiles_per_window,eval_ebn0_db,stage1_bits,stage2_bits,keep_ratio,"
          "interleaver_name,bits_per_symbol,normalize_extrinsic,generate_random_bits,"
          "normalize_known_prefix_tail,llr_bits,quant_clip_ratio,siso_active_list,"
          "alpha_low_grid,alpha_high_grid,beta_low_grid,beta_high_grid,gamma_alpha_grid,gamma_beta_grid,"
          "decoder_name,"
          "alpha_low,alpha_high,gamma_alpha,beta_low,beta_high,gamma_beta,"
          "early_stop_beta_start,early_stop_beta_step,"
          "early_stop_action_hard_llr_mag,"
          "chase_L,chase_n_test,chase_topk_keep,chase_group_minima_bits,bitgen_seed,channel_seed,ebn0_db,"
          "alpha_list,beta_list,early_stop_beta_list,"
          "early_stop_condition_mode,early_stop_action_mode,early_stop_bind_group_size,"
          "early_stop_condition_mode_list_cfg,early_stop_action_mode_list_cfg,early_stop_bind_group_size_list_cfg,"
          "early_stop_cond_v1_require_bch,early_stop_cond_v1_require_overall,"
          "early_stop_v2_llr_abs_threshold,early_stop_v2_max_unreliable_bits,"
          "early_stop_cond_v2_include_overall,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "early_stop_list,early_stop_row_list,"
          "unscheduled_count_list,unscheduled_list,unscheduled_need_list\n";
}

std::string join_vec(const std::vector<float>& values, char sep, int precision) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss << std::setprecision(precision);
  for (size_t i = 0; i < values.size(); ++i) {
    oss << values[i];
    if (i + 1 < values.size()) {
      oss << sep;
    }
  }
  return oss.str();
}

std::string join_vec(const std::vector<double>& values, char sep, int precision) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss << std::setprecision(precision);
  for (size_t i = 0; i < values.size(); ++i) {
    oss << values[i];
    if (i + 1 < values.size()) {
      oss << sep;
    }
  }
  return oss.str();
}

void write_csv_row(std::ostream& csv,
                   const std::string& timestamp,
                   const std::string& run_id,
                   const std::string& stage_tag,
                   std::size_t num_bits,
                   const SweepScenario& scenario,
                   const newcode::PipelineResult& result,
                   CsvFormat format,
                   const ExtendedCsvConfigSnapshot* snapshot) {
  const auto prev_prec = csv.precision();
  if (format == CsvFormat::Basic) {
    const ExtendedCsvConfigSnapshot empty_snapshot{};
    const ExtendedCsvConfigSnapshot& cfg =
        snapshot ? *snapshot : empty_snapshot;
    const double es_mean = mean(result.tile_early_stop_pct);
    const double es_row_mean = mean(result.tile_row_early_stop_pct);
    const double unsched_mean = mean(result.tile_unscheduled_pct);
    const double unsched_need_mean = mean(result.tile_unscheduled_among_need_pct);
    csv << timestamp << ","
        << run_id << ","
        << scenario.name << ","
        << scenario.decoder_name << ","
        << result.pre_fec.ber << ","
        << result.pre_fec.errors << ","
        << result.pre_fec.total << ","
        << std::setprecision(10) << result.post_fec.ber << ","
        << result.post_fec.errors << ","
        << result.post_fec.total << ","
        << scenario.early_stop_action_hard_llr_mag << ","
        << scenario.chase_L << ","
        << scenario.chase_n_test << ","
        << scenario.chase_topk_keep << ","
        << scenario.chase_group_minima_bits << ","
        << scenario.mux_group_g << ","
        << scenario.mux_scheduling_mode << ","
        << scenario.mux_early_stop_priority_rule << ","
        << scenario.mux_bypass_scheme << ","
        << scenario.bitgen_seed << ","
        << scenario.channel_seed << ","
        << scenario.ebn0_db << ","
        << '"' << join_vec(scenario.alpha_list, '|', 6) << "\","
        << '"' << join_vec(scenario.beta_list, '|', 6) << "\","
        << '"' << join_vec(scenario.early_stop_action_sign_beta_list, '|', 6) << "\","
        << scenario.early_stop_condition_mode << ","
        << scenario.early_stop_action_mode << ","
        << scenario.early_stop_bind_group_size << ","
        << (scenario.early_stop_cond_v1_require_bch ? 1 : 0) << ","
        << (scenario.early_stop_cond_v1_require_overall ? 1 : 0) << ","
        << scenario.early_stop_v2_llr_abs_threshold << ","
        << scenario.early_stop_v2_max_unreliable_bits << ","
        << (scenario.early_stop_cond_v2_include_overall ? 1 : 0) << ",";
    csv.precision(prev_prec);
    if (std::isnan(es_mean)) {
      csv << ",";
    } else {
      csv << std::setprecision(3) << es_mean << ",";
      csv.precision(prev_prec);
    }
    if (std::isnan(es_row_mean)) {
      csv << ",";
    } else {
      csv << std::setprecision(3) << es_row_mean << ",";
      csv.precision(prev_prec);
    }
    csv << '"' << join_vec(result.tile_early_stop_pct, '|', 1) << "\","
        << '"' << join_vec(result.tile_row_early_stop_pct, '|', 1) << "\",";
    if (std::isnan(unsched_mean)) {
      csv << ",";
    } else {
      csv << std::setprecision(3) << unsched_mean << ",";
      csv.precision(prev_prec);
    }
    if (std::isnan(unsched_need_mean)) {
      csv << ",";
    } else {
      csv << std::setprecision(3) << unsched_need_mean << ",";
      csv.precision(prev_prec);
    }
    csv << '"' << join_vec_size_t(result.tile_unscheduled_count, '|') << "\","
        << '"' << join_vec(result.tile_unscheduled_pct, '|', 1) << "\","
        << '"' << join_vec(result.tile_unscheduled_among_need_pct, '|', 1) << "\","
        << '"' << cfg.interleaver_name << '"' << ","
        << cfg.bits_per_symbol << ","
        << (cfg.generate_random_bits ? 1 : 0) << ","
        << cfg.bitgen_seed << ","
        << cfg.bitgen_seed_count << ","
        << '"' << cfg.bitgen_seed_candidates << '"' << ","
        << cfg.ebn0_start << ","
        << cfg.ebn0_end << ","
        << cfg.ebn0_points << ","
        << '"' << cfg.ebn0_candidates << '"' << ","
        << cfg.channel_seed << ","
        << cfg.channel_seed_count << ","
        << '"' << cfg.channel_seed_candidates << '"' << ","
        << cfg.llr_bits << ","
        << cfg.quant_clip_ratio << ","
        << (cfg.normalize_extrinsic ? 1 : 0) << ","
        << (cfg.normalize_known_prefix_tail ? 1 : 0) << ","
        << '"' << cfg.decoder_name_candidates << '"' << ","
        << '"' << cfg.chase_l_candidates << '"' << ","
        << cfg.chase_n_test << ","
        << '"' << cfg.chase_n_test_candidates << '"' << ","
        << cfg.chase_topk_keep << ","
        << '"' << cfg.chase_topk_keep_candidates << '"' << ","
        << cfg.chase_group_minima_bits << ","
        << '"' << cfg.chase_group_minima_bits_candidates << '"' << ","
        << '"' << cfg.alpha_start_candidates << '"' << ","
        << '"' << cfg.alpha_step_candidates << '"' << ","
        << '"' << cfg.beta_start_candidates << '"' << ","
        << '"' << cfg.beta_step_candidates << '"' << ","
        << '"' << cfg.explicit_patterns << '"' << ","
        << '"' << join_vec_int(cfg.siso_active_list, '|') << '"' << ","
        << cfg.mux_group_g << ","
        << cfg.mux_scheduling_mode << ","
        << '"' << cfg.mux_scheduling_mode_candidates << '"' << ","
        << cfg.mux_early_stop_priority_rule << ","
        << '"' << cfg.mux_early_stop_priority_rule_candidates << '"' << ","
        << (cfg.mux_enable_reconfig ? 1 : 0) << ","
        << cfg.mux_bypass_scheme << ","
        << (cfg.enable_early_stop ? 1 : 0) << ","
        << cfg.early_stop_condition_mode << ","
        << '"' << cfg.early_stop_condition_candidates << '"' << ","
        << '"' << cfg.early_stop_condition_mode_list << '"' << ","
        << cfg.early_stop_action_mode << ","
        << '"' << cfg.early_stop_action_mode_list << '"' << ","
        << cfg.early_stop_bind_group_size << ","
        << '"' << cfg.early_stop_bind_group_size_list << '"' << ","
        << '"' << cfg.early_stop_action_candidates << '"' << ","
        << (cfg.early_stop_cond_v1_require_bch ? 1 : 0) << ","
        << '"' << cfg.early_stop_cond_v1_require_bch_candidates << '"' << ","
        << (cfg.early_stop_cond_v1_require_overall ? 1 : 0) << ","
        << '"' << cfg.early_stop_cond_v1_require_overall_candidates << '"' << ","
        << cfg.early_stop_v2_llr_abs_threshold << ","
        << '"' << cfg.early_stop_v2_llr_abs_threshold_candidates << '"' << ","
        << cfg.early_stop_v2_max_unreliable_bits << ","
        << '"' << cfg.early_stop_v2_max_unreliable_bits_candidates << '"' << ","
        << (cfg.early_stop_cond_v2_include_overall ? 1 : 0) << ","
        << cfg.early_stop_action_sign_beta_fill << ","
        << '"' << cfg.early_stop_action_beta_start_candidates << '"' << ","
        << '"' << cfg.early_stop_action_beta_step_candidates << '"' << ","
        << cfg.early_stop_action_residual_divisor << ","
        << cfg.early_stop_action_hard_llr_mag << ","
        << '"' << cfg.early_stop_action_hard_llr_mag_candidates << '"' << ","
        << (cfg.quiet_pipeline ? 1 : 0) << ","
        << (cfg.quiet_logs ? 1 : 0) << ","
        << (cfg.trace_enable ? 1 : 0) << ","
        << cfg.trace_row << ","
        << cfg.trace_col << ","
        << (cfg.trace_log_read ? 1 : 0) << ","
        << (cfg.trace_log_write ? 1 : 0) << ","
        << (cfg.trace_log_mismatch ? 1 : 0) << "\n";
  } else {
    const ExtendedCsvConfigSnapshot empty_snapshot{};
    const ExtendedCsvConfigSnapshot& cfg =
        snapshot ? *snapshot : empty_snapshot;
    csv << timestamp << ","
        << run_id << ","
        << stage_tag << ","
        << num_bits << ","
        << '"' << scenario.name << "\","
        << cfg.tiles_per_window << ","
        << cfg.eval_ebn0_db << ","
        << cfg.stage1_bits << ","
        << cfg.stage2_bits << ","
        << cfg.keep_ratio << ","
        << '"' << cfg.interleaver_name << '"' << ","
        << cfg.bits_per_symbol << ","
        << (cfg.normalize_extrinsic ? 1 : 0) << ","
        << (cfg.generate_random_bits ? 1 : 0) << ","
        << (cfg.normalize_known_prefix_tail ? 1 : 0) << ","
        << cfg.llr_bits << ","
        << cfg.quant_clip_ratio << ","
        << '"' << join_vec_int(cfg.siso_active_list, '|') << '"' << ","
        << '"' << join_vec(cfg.alpha_low_grid, '|', 6) << '"' << ","
        << '"' << join_vec(cfg.alpha_high_grid, '|', 6) << '"' << ","
        << '"' << join_vec(cfg.beta_low_grid, '|', 6) << '"' << ","
        << '"' << join_vec(cfg.beta_high_grid, '|', 6) << '"' << ","
        << '"' << join_vec(cfg.gamma_alpha_grid, '|', 6) << '"' << ","
        << '"' << join_vec(cfg.gamma_beta_grid, '|', 6) << '"' << ","
        << scenario.decoder_name << ","
        << scenario.alpha_low << ","
        << scenario.alpha_high << ","
        << scenario.gamma_alpha << ","
        << scenario.beta_low << ","
        << scenario.beta_high << ","
        << scenario.gamma_beta << ","
        << scenario.early_stop_beta_start << ","
        << scenario.early_stop_beta_step << ","
        << scenario.early_stop_action_hard_llr_mag << ","
        << scenario.chase_L << ","
        << scenario.chase_n_test << ","
        << scenario.chase_topk_keep << ","
        << scenario.chase_group_minima_bits << ","
        << scenario.bitgen_seed << ","
        << scenario.channel_seed << ","
        << scenario.ebn0_db << ","
        << '"' << join_vec(scenario.alpha_list, '|', 6) << "\","
        << '"' << join_vec(scenario.beta_list, '|', 6) << "\","
        << '"' << join_vec(scenario.early_stop_action_sign_beta_list, '|', 6) << "\","
        << scenario.early_stop_condition_mode << ","
        << scenario.early_stop_action_mode << ","
        << scenario.early_stop_bind_group_size << ","
        << '"' << cfg.early_stop_condition_mode_list << '"' << ","
        << '"' << cfg.early_stop_action_mode_list << '"' << ","
        << '"' << cfg.early_stop_bind_group_size_list << '"' << ","
        << (scenario.early_stop_cond_v1_require_bch ? 1 : 0) << ","
        << (scenario.early_stop_cond_v1_require_overall ? 1 : 0) << ","
        << scenario.early_stop_v2_llr_abs_threshold << ","
        << scenario.early_stop_v2_max_unreliable_bits << ","
        << (scenario.early_stop_cond_v2_include_overall ? 1 : 0) << ","
        << result.pre_fec.ber << ","
        << result.pre_fec.errors << ","
        << result.pre_fec.total << ","
        << std::setprecision(10) << result.post_fec.ber << ","
        << result.post_fec.errors << ","
        << result.post_fec.total << ","
        << '"' << join_vec(result.tile_early_stop_pct, '|', 1) << "\","
        << '"' << join_vec(result.tile_row_early_stop_pct, '|', 1) << "\","
        << '"' << join_vec_size_t(result.tile_unscheduled_count, '|') << "\","
        << '"' << join_vec(result.tile_unscheduled_pct, '|', 1) << "\","
        << '"' << join_vec(result.tile_unscheduled_among_need_pct, '|', 1) << "\","
        << '"' << cfg.interleaver_name << '"' << ","
        << cfg.bits_per_symbol << ","
        << (cfg.generate_random_bits ? 1 : 0) << ","
        << cfg.bitgen_seed << ","
        << cfg.bitgen_seed_count << ","
        << '"' << cfg.bitgen_seed_candidates << '"' << ","
        << cfg.ebn0_start << ","
        << cfg.ebn0_end << ","
        << cfg.ebn0_points << ","
        << '"' << cfg.ebn0_candidates << '"' << ","
        << cfg.channel_seed << ","
        << cfg.channel_seed_count << ","
        << '"' << cfg.channel_seed_candidates << '"' << ","
        << cfg.llr_bits << ","
        << cfg.quant_clip_ratio << ","
        << (cfg.normalize_extrinsic ? 1 : 0) << ","
        << (cfg.normalize_known_prefix_tail ? 1 : 0) << ","
        << '"' << cfg.decoder_name_candidates << '"' << ","
        << '"' << cfg.chase_l_candidates << '"' << ","
        << cfg.chase_n_test << ","
        << '"' << cfg.chase_n_test_candidates << '"' << ","
        << cfg.chase_topk_keep << ","
        << '"' << cfg.chase_topk_keep_candidates << '"' << ","
        << cfg.chase_group_minima_bits << ","
        << '"' << cfg.chase_group_minima_bits_candidates << '"' << ","
        << '"' << cfg.alpha_start_candidates << '"' << ","
        << '"' << cfg.alpha_step_candidates << '"' << ","
        << '"' << cfg.beta_start_candidates << '"' << ","
        << '"' << cfg.beta_step_candidates << '"' << ","
        << '"' << cfg.explicit_patterns << '"' << ","
        << '"' << join_vec_int(cfg.siso_active_list, '|') << '"' << ","
        << cfg.mux_group_g << ","
        << (cfg.mux_enable_reconfig ? 1 : 0) << ","
        << cfg.mux_bypass_scheme << ","
        << (cfg.enable_early_stop ? 1 : 0) << ","
        << cfg.early_stop_condition_mode << ","
        << '"' << cfg.early_stop_condition_candidates << '"' << ","
        << cfg.early_stop_action_mode << ","
        << '"' << cfg.early_stop_action_candidates << '"' << ","
        << (cfg.early_stop_cond_v1_require_bch ? 1 : 0) << ","
        << '"' << cfg.early_stop_cond_v1_require_bch_candidates << '"' << ","
        << (cfg.early_stop_cond_v1_require_overall ? 1 : 0) << ","
        << '"' << cfg.early_stop_cond_v1_require_overall_candidates << '"' << ","
        << cfg.early_stop_v2_llr_abs_threshold << ","
        << '"' << cfg.early_stop_v2_llr_abs_threshold_candidates << '"' << ","
        << cfg.early_stop_v2_max_unreliable_bits << ","
        << '"' << cfg.early_stop_v2_max_unreliable_bits_candidates << '"' << ","
        << (cfg.early_stop_cond_v2_include_overall ? 1 : 0) << ","
        << cfg.early_stop_action_sign_beta_fill << ","
        << '"' << cfg.early_stop_action_beta_start_candidates << '"' << ","
        << '"' << cfg.early_stop_action_beta_step_candidates << '"' << ","
        << cfg.early_stop_action_residual_divisor << ","
        << cfg.early_stop_action_hard_llr_mag << ","
        << '"' << cfg.early_stop_action_hard_llr_mag_candidates << '"' << ","
        << (cfg.quiet_pipeline ? 1 : 0) << ","
        << (cfg.quiet_logs ? 1 : 0) << ","
        << (cfg.trace_enable ? 1 : 0) << ","
        << cfg.trace_row << ","
        << cfg.trace_col << ","
        << (cfg.trace_log_read ? 1 : 0) << ","
        << (cfg.trace_log_write ? 1 : 0) << ","
        << (cfg.trace_log_mismatch ? 1 : 0) << "\n";
    csv.precision(prev_prec);
  }
}

}  // namespace detail
}  // namespace ofec_sweep
