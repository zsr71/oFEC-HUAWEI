#include "newcode/ofec_single_runner.hpp"

namespace ofec_single {

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

}  // namespace ofec_single
