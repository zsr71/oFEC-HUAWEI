#include "tpc_sweep_detail.hpp"

#include <iomanip>
#include <sstream>

namespace tpc_sweep {
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
  fout << "timestamp,run_id,scenario,bitgen_seed,channel_seed,ebn0_db,"
          "bits_per_symbol,max_iters,num_blocks,alpha_schedule,beta_schedule,"
          "pre_ber,pre_errs,pre_total,post_ber,post_errs,post_total,"
          "avg_iters,early_stop_blocks,blocks\n";
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

void write_csv_row(std::ostream& csv,
                   const std::string& timestamp,
                   const std::string& run_id,
                   const SweepScenario& scenario,
                   const newcode::tpc::TpcPipelineResult& result,
                   const SweepParameterConfig& config) {
  const auto prev_prec = csv.precision();
  csv << timestamp << ","
      << run_id << ","
      << scenario.name << ","
      << scenario.bitgen_seed << ","
      << scenario.channel_seed << ","
      << scenario.ebn0_db << ","
      << config.bits_per_symbol << ","
      << config.max_iters << ","
      << config.num_blocks << ","
      << '"' << join_vec(config.alpha_schedule, '|', 3) << "\","
      << '"' << join_vec(config.beta_schedule, '|', 3) << "\","
      << result.pre_fec.ber << ","
      << result.pre_fec.errors << ","
      << result.pre_fec.total << ","
      << std::setprecision(10) << result.post_fec.ber << ","
      << result.post_fec.errors << ","
      << result.post_fec.total << ","
      << std::setprecision(3) << result.avg_iters << ","
      << result.early_stop_blocks << ","
      << result.blocks << "\n";
  csv.precision(prev_prec);
}

}  // namespace detail
}  // namespace tpc_sweep
