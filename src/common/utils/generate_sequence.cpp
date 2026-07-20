#include "ofec_sweep_detail.hpp"

namespace ofec_sweep {
namespace detail {

std::vector<float> generate_sequence(float start, float step, std::size_t length) {
  std::vector<float> seq(length, start);
  for (std::size_t i = 0; i < length; ++i) {
    seq[i] = start + step * static_cast<float>(i);
  }
  return seq;
}

}  // namespace detail
}  // namespace ofec_sweep
