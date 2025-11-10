#include "newcode/ofec_single_runner.hpp"

namespace ofec_single {

DualWriter::DualWriter(std::ofstream& file)
    : console_(&std::cout), file_(&file) {}

DualWriter& DualWriter::operator<<(std::ostream& (*manip)(std::ostream&)) {
  if (console_) {
    manip(*console_);
  }
  if (file_ && file_->is_open()) {
    manip(*file_);
  }
  return *this;
}

DualWriter& DualWriter::operator<<(std::ios_base& (*manip)(std::ios_base&)) {
  if (console_) {
    manip(*console_);
  }
  if (file_ && file_->is_open()) {
    manip(*file_);
  }
  return *this;
}

int run_ofec_single(const Config& config) {
  const std::filesystem::path data_dir = "data";
  std::string log_path;
  std::ofstream log_file = detail::prepare_log_file(data_dir, log_path);
  DualWriter log(log_file);

  std::optional<newcode::Params> params = detail::build_params(config, log);
  if (!params.has_value()) {
    return 2;
  }

  detail::log_run_overview(config, *params, log);
  newcode::PipelineConfig pipeline_cfg = detail::build_pipeline_config(config);
  const newcode::PipelineResult result =
      newcode::run_pipeline(*params, pipeline_cfg, config.label, config.ebn0_db);
  detail::log_pipeline_results(result, log);
  log << "[INFO] log saved at " << log_path << "\n";
  return 0;
}

}  // namespace ofec_single
