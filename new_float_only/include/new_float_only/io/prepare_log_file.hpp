#pragma once

#include <filesystem>
#include <fstream>
#include <string>

namespace io {

/**
 * 在 data 目录下创建单次运行日志文件。
 * log_path 会返回最终文件路径，返回值为已打开的输出文件流。
 */
std::ofstream prepare_log_file(const std::filesystem::path& data_dir,
                               std::string& log_path);

}  // namespace io
