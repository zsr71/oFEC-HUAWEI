#include "new_float_only/io/ensure_dir.hpp"
#include "new_float_only/io/prepare_log_file.hpp"
#include "new_float_only/utils/now_stamp.hpp"

namespace io {

/**
 * 创建当前单次运行使用的日志文件。
 * 关键语句：先确保目录存在，再按时间戳拼接文件名，最后以 append 模式打开。
 */
std::ofstream prepare_log_file(const std::filesystem::path& data_dir,
                               std::string& log_path) {
  ensure_dir(data_dir);
  log_path = (data_dir / ("run_" + utils::now_stamp() + "_single.log")).string();
  std::ofstream file(log_path, std::ios::out | std::ios::app);
  return file;
}

}  // namespace io
