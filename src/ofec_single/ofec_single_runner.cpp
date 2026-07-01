#include "newcode/ofec_single_runner.hpp"
#include "newcode/io/dualwriter.hpp"
#include "newcode/io/prepare_log_file.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

namespace {

std::string make_run_id() {
  const auto now = std::chrono::system_clock::now();
  const std::time_t now_time = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &now_time);
#else
  localtime_r(&now_time, &tm);
#endif
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
  return oss.str();
}

void dump_tile_early_stop_samples_csv(
    const std::vector<newcode::TileEarlyStopSample>& samples,
    const std::filesystem::path& output_path) {
  std::filesystem::path parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "invocation,tile_index,rows_total,rows_passed,"
         "rows_hard_finish,rows_need_siso_before_mux,rows_unscheduled\n";
  for (const auto& sample : samples) {
    out << sample.invocation << ','
        << sample.tile_index << ','
        << sample.rows_total << ','
        << sample.rows_passed << ','
        << sample.rows_hard_finish << ','
        << sample.rows_need_siso_before_mux << ','
        << sample.rows_unscheduled << '\n';
  }
}

void dump_tile_early_stop_group_bind_debug_samples_csv(
    const std::vector<newcode::TileEarlyStopGroupBindDebugSample>& samples,
    const std::filesystem::path& output_path,
    const std::string& run_id,
    const std::string& label) {
  std::filesystem::path parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "run_id,label,invocation,tile_index,condition_mode,bind_group_size,"
         "raw_early_stop_flags,bound_early_stop_flags\n";
  for (const auto& sample : samples) {
    out << run_id << ','
        << label << ','
        << sample.invocation << ','
        << sample.tile_index << ','
        << sample.condition_mode << ','
        << sample.bind_group_size << ','
        << sample.raw_early_stop_flags << ','
        << sample.bound_early_stop_flags << '\n';
  }
}

}  // namespace

namespace ofec_single {

int run_ofec_single(const Config& config) {
  const std::filesystem::path data_dir = "data";
  const std::string run_id = make_run_id();
  std::string log_path;
  std::ofstream log_file = io::prepare_log_file(data_dir, log_path);
  io::DualWriter log(log_file);

  std::optional<newcode::Params> params = detail::build_params(config, log);
  if (!params.has_value()) {
    return 2;
  }

  detail::log_run_overview(config, *params, log);
  newcode::PipelineConfig pipeline_cfg = detail::build_pipeline_config(config);
  const newcode::PipelineResult result =
      newcode::run_pipeline(*params, pipeline_cfg, config.label, config.ebn0_db);
  if (config.dump_tile_early_stop_samples) {
    std::filesystem::path output_path = config.tile_early_stop_samples_output_path;
    if (output_path.empty()) {
      output_path = std::filesystem::path("data/early_stop_hist") /
                    (config.label + "_tile_early_stop_samples.csv");
    }
    dump_tile_early_stop_samples_csv(result.tile_early_stop_samples, output_path);
    log << "[INFO] tile early-stop samples saved to " << output_path.string() << "\n";
  }
  if (config.dump_tile_early_stop_group_bind_debug_samples) {
    std::filesystem::path output_path =
        config.tile_early_stop_group_bind_debug_samples_output_path;
    if (output_path.empty()) {
      output_path = std::filesystem::path("data/early_stop_debug") /
                    (config.label + "_group_bind_debug.csv");
    }
    dump_tile_early_stop_group_bind_debug_samples_csv(
        result.tile_early_stop_group_bind_debug_samples,
        output_path,
        run_id,
        config.label);
    log << "[INFO] tile early-stop group-bind debug samples saved to "
        << output_path.string() << "\n";
  }
  detail::log_pipeline_results(result, log);
  log << "[INFO] log saved at " << log_path << "\n";
  return 0;
}

}  // namespace ofec_single
