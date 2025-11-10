#include "ofec_sweep_detail.hpp"

#include <iomanip>
#include <limits>
#include <sstream>
#include <system_error>

namespace ofec_sweep {
namespace detail {

DualOut::DualOut(std::ostream& console, const std::string& filepath)
    : console_(console), file_(filepath, std::ios::out | std::ios::app) {}

DualOut& DualOut::operator<<(std::ostream& (*pf)(std::ostream&)) {
  pf(console_);
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
          "early_stop_mean_pct,early_stop_list\n";
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

}  // namespace detail
}  // namespace ofec_sweep

