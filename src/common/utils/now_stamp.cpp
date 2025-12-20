#include "newcode/utils/now_stamp.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>

namespace utils {

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

}  // namespace utils
