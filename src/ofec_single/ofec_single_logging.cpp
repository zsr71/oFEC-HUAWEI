#include "newcode/ofec_single_runner.hpp"

#include <chrono>
#include <iomanip>
#include <sstream>
#include <system_error>

namespace ofec_single {
namespace {

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

void ensure_dir(const std::filesystem::path& p) {
  std::error_code ec;
  std::filesystem::create_directories(p, ec);
}

}  // namespace

namespace detail {

std::ofstream prepare_log_file(const std::filesystem::path& data_dir,
                               std::string& log_path) {
  ensure_dir(data_dir);
  log_path = (data_dir / ("run_" + now_stamp() + "_single.log")).string();
  std::ofstream file(log_path, std::ios::out | std::ios::app);
  return file;
}

}  // namespace detail
}  // namespace ofec_single
