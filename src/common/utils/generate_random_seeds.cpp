#include "ofec_sweep_detail.hpp"

#include <limits>
#include <random>
#include <unordered_set>

namespace ofec_sweep {
namespace detail {

std::vector<int> generate_random_seeds(int count) {
  std::vector<int> seeds;
  if (count <= 0) {
    return seeds;
  }

  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<int> dist(1, std::numeric_limits<int>::max());
  std::unordered_set<int> seen;
  seen.reserve(static_cast<std::size_t>(count));
  seeds.reserve(static_cast<std::size_t>(count));

  while (seeds.size() < static_cast<std::size_t>(count)) {
    int candidate = dist(rng);
    if (seen.insert(candidate).second) {
      seeds.push_back(candidate);
    }
  }
  return seeds;
}

}  // namespace detail
}  // namespace ofec_sweep
