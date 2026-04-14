#include "newcode/ofec_single_runner.hpp"
#include "newcode/io/dualwriter.hpp"
#include "newcode/io/prepare_log_file.hpp"

namespace {

void dump_tile_early_stop_samples_csv(
    const std::vector<newcode::TileEarlyStopSample>& samples,
    const std::filesystem::path& output_path) {
  std::filesystem::path parent = output_path.parent_path();
  if (!parent.empty()) {
    std::filesystem::create_directories(parent);
  }
  std::ofstream out(output_path);
  out << "invocation,tile_index,rows_total,rows_passed\n";
  for (const auto& sample : samples) {
    out << sample.invocation << ','
        << sample.tile_index << ','
        << sample.rows_total << ','
        << sample.rows_passed << '\n';
  }
}

}  // namespace

namespace ofec_single {

int run_ofec_single(const Config& config) {
  const std::filesystem::path data_dir = "data";
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
  detail::log_pipeline_results(result, log);
  log << "[INFO] log saved at " << log_path << "\n";
  return 0;
}

}  // namespace ofec_single
