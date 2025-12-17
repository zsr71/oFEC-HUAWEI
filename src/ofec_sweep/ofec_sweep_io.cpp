#include "ofec_sweep_detail.hpp"

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <system_error>

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

std::string now_stamp() {
  using clock = std::chrono::system_clock;
  const auto t = clock::to_time_t(clock::now());
  std::tm tm{};
#ifdef _WIN32
  localtime_s(&tm, &t);
#else
  localtime_r(&t, &tm);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d-%H%M%S");
  return oss.str();
}

void ensure_dir(const std::filesystem::path& path) {
  std::error_code ec;
  std::filesystem::create_directories(path, ec);
}

void ensure_csv_header(const std::string& csv_path) {
  std::ifstream fin(csv_path);
  if (fin.good() && fin.peek() != std::ifstream::traits_type::eof()) {
    return;
  }

  std::ofstream fout(csv_path, std::ios::out | std::ios::app);
  fout << "timestamp,run_id,scenario,alpha_start,alpha_step,beta_start,beta_step,"
          "chase_L,chase_n_test,bitgen_seed,channel_seed,ebn0_db,ALPHA_LIST,beta_list,"
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
          "alpha_low,alpha_high,gamma_alpha,beta_low,beta_high,gamma_beta,"
          "chase_L,chase_n_test,bitgen_seed,channel_seed,ebn0_db,"
          "alpha_list,beta_list,pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
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

double mean(const std::vector<double>& values) {
  if (values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  long double sum = 0.0L;
  for (double v : values) {
    sum += v;
  }
  return static_cast<double>(sum / values.size());
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
        << scenario.alpha_start << ","
        << alpha_step << ","
        << scenario.beta_start << ","
        << beta_step << ","
        << scenario.chase_L << ","
        << scenario.chase_n_test << ","
        << scenario.bitgen_seed << ","
        << scenario.channel_seed << ","
        << scenario.ebn0_db << ","
        << '"' << join_vec(scenario.alpha_list, '|', 6) << "\","
        << '"' << join_vec(scenario.beta_list, '|', 6) << "\","
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
        << scenario.alpha_low << ","
        << scenario.alpha_high << ","
        << scenario.gamma_alpha << ","
        << scenario.beta_low << ","
        << scenario.beta_high << ","
        << scenario.gamma_beta << ","
        << scenario.chase_L << ","
        << scenario.chase_n_test << ","
        << scenario.bitgen_seed << ","
        << scenario.channel_seed << ","
        << scenario.ebn0_db << ","
        << '"' << join_vec(scenario.alpha_list, '|', 6) << "\","
        << '"' << join_vec(scenario.beta_list, '|', 6) << "\","
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
