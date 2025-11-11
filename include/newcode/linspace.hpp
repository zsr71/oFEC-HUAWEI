#pragma once

#include <cstddef>
#include <initializer_list>
#include <vector>

namespace newcode {

struct LinspaceSpec {
  float start;
  float end;
  std::size_t count;
};

std::vector<float> linspace(float start, float end, std::size_t count);


}  // namespace newcode

