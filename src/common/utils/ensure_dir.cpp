#include "newcode/ensure_dir.hpp"

#include <system_error>

namespace newcode {

void ensure_dir(const std::filesystem::path& p) {
  std::error_code ec;
  std::filesystem::create_directories(p, ec);
}

}  // namespace newcode
