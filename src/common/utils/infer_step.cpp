#include "ofec_sweep_detail.hpp"

namespace ofec_sweep {
namespace detail {

float infer_step(const std::vector<float>& values) {
  return values.size() >= 2 ? values[1] - values[0] : 0.0f;
}

}  // namespace detail
}  // namespace ofec_sweep
