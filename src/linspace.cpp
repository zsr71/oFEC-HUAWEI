#include "newcode/linspace.hpp"

namespace newcode {

std::vector<float> linspace(float start, float end, std::size_t count) {
  std::vector<float> values;
  if (count == 0) {
    return values;
  }
  if (count == 1) {
    values.push_back(start);
    return values;
  }
  values.reserve(count);
  const float step = (end - start) / static_cast<float>(count - 1);
  for (std::size_t i = 0; i < count; ++i) {
    values.push_back(start + step * static_cast<float>(i));
  }
  return values;
}

}  // namespace newcode

