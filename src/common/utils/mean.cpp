#include "ofec_sweep_detail.hpp"

#include <limits>

namespace ofec_sweep {
namespace detail {

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
