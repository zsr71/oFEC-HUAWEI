#include "ofec_sweep_detail.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace ofec_sweep {
namespace detail {

DualOut::DualOut(std::ostream& console, const std::string& filepath, bool mirror_console)
    : console_(mirror_console ? &console : nullptr),
      file_(filepath, std::ios::out | std::ios::app) {}

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
  fout << "timestamp,run_id,scenario,decoder_name,alpha_start,alpha_step,beta_start,beta_step,"
          "early_stop_beta_start,early_stop_beta_step,"
          "early_stop_action_hard_llr_mag,"
          "chase_L,chase_n_test,chase_topk_keep,chase_group_minima_bits,bitgen_seed,channel_seed,ebn0_db,ALPHA_LIST,beta_list,"
          "early_stop_beta_list,early_stop_condition_mode,early_stop_action_mode,"
          "early_stop_cond_v1_require_bch,early_stop_cond_v1_require_overall,"
          "early_stop_v2_llr_abs_threshold,early_stop_v2_max_unreliable_bits,"
          "early_stop_cond_v2_include_overall,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "early_stop_mean_pct,early_stop_row_mean_pct,early_stop_list,early_stop_row_list\n";
}

void ensure_csv_header_v2(const std::string& csv_path) {
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,stage,num_bits,label,"
          "decoder_name,"
          "alpha_low,alpha_high,gamma_alpha,beta_low,beta_high,gamma_beta,"
          "early_stop_beta_start,early_stop_beta_step,"
          "early_stop_action_hard_llr_mag,"
          "chase_L,chase_n_test,chase_topk_keep,chase_group_minima_bits,bitgen_seed,channel_seed,ebn0_db,"
          "alpha_list,beta_list,early_stop_beta_list,"
          "early_stop_condition_mode,early_stop_action_mode,"
          "early_stop_cond_v1_require_bch,early_stop_cond_v1_require_overall,"
          "early_stop_v2_llr_abs_threshold,early_stop_v2_max_unreliable_bits,"
          "early_stop_cond_v2_include_overall,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "early_stop_list,early_stop_row_list\n";
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
                   CsvFormat format) {
  const auto prev_prec = csv.precision();
  if (format == CsvFormat::Basic) {
    const float alpha_step = infer_step(scenario.alpha_list);
    const float beta_step = infer_step(scenario.beta_list);
    const double es_mean = mean(result.tile_early_stop_pct);
    const double es_row_mean = mean(result.tile_row_early_stop_pct);
    csv << timestamp << ","
        << run_id << ","
        << scenario.name << ","
        << scenario.decoder_name << ","
        << scenario.alpha_start << ","
        << alpha_step << ","
        << scenario.beta_start << ","
        << beta_step << ","
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
        << result.post_fec.total << ",";
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
        << '"' << join_vec(result.tile_row_early_stop_pct, '|', 1) << "\"\n";
  } else {
    csv << timestamp << ","
        << run_id << ","
        << stage_tag << ","
        << num_bits << ","
        << '"' << scenario.name << "\","
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
        << '"' << join_vec(result.tile_row_early_stop_pct, '|', 1) << "\"\n";
    csv.precision(prev_prec);
  }
}

}  // namespace detail
}  // namespace ofec_sweep
