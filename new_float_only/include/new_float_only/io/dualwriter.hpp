
#pragma once
#include <filesystem>
#include <fstream>
#include <iostream>


namespace io{

class DualWriter {
 public:
  explicit DualWriter(std::ofstream& file);

  template <typename T>
  DualWriter& operator<<(const T& value) {
    if (console_) {
      *console_ << value;
    }
    if (file_ && file_->is_open()) {
      *file_ << value;
    }
    return *this;
  }

  DualWriter& operator<<(std::ostream& (*manip)(std::ostream&));
  DualWriter& operator<<(std::ios_base& (*manip)(std::ios_base&));

 private:
  std::ostream* console_;
  std::ofstream* file_;
};

}