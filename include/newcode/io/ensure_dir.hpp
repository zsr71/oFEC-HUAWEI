#pragma once

#include <filesystem>

namespace io {

void ensure_dir(const std::filesystem::path& path);

}  // namespace io