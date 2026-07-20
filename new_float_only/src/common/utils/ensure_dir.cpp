#include "new_float_only/io/ensure_dir.hpp"

#include <system_error>

namespace io {

void ensure_dir(const std::filesystem::path& p) {
  std::error_code ec;
  std::filesystem::create_directories(p, ec);
}

}  // namespace new_float_only
