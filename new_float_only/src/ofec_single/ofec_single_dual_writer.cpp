#include "new_float_only/io/dualwriter.hpp"

namespace io {

/**
 * DualWriter 同时把日志写到控制台和文件。
 * 构造时绑定外部已经打开的文件流，后续所有输出都复用这两个目标。
 */
DualWriter::DualWriter(std::ofstream& file)
    : console_(&std::cout), file_(&file) {}

DualWriter& DualWriter::operator<<(std::ostream& (*manip)(std::ostream&)) {
  if (console_) {
    manip(*console_);
  }
  if (file_ && file_->is_open()) {
    manip(*file_);
  }
  return *this;
}

DualWriter& DualWriter::operator<<(std::ios_base& (*manip)(std::ios_base&)) {
  if (console_) {
    manip(*console_);
  }
  if (file_ && file_->is_open()) {
    manip(*file_);
  }
  return *this;
}

}  // namespace io
