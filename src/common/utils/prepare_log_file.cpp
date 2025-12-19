#include "newcode/ofec_single_runner.hpp"
#include "newcode/ensure_dir.hpp"
#include "newcode/now_stamp.hpp"

namespace io {


std::ofstream prepare_log_file(const std::filesystem::path& data_dir,
                               std::string& log_path) {
  newcode::ensure_dir(data_dir);
  log_path = (data_dir / ("run_" + newcode::now_stamp() + "_single.log")).string();
  std::ofstream file(log_path, std::ios::out | std::ios::app);
  return file;


}  // namespace io
}