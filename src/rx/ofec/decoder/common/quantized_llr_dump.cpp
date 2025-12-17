#include "newcode/quantized_llr_dump.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <system_error>
#include <vector>

namespace newcode {
namespace {

template <typename T>
std::vector<T> flatten_matrix(const Matrix<T>& matrix) {
  std::vector<T> out;
  out.reserve(matrix.rows() * matrix.cols());
  for (size_t r = 0; r < matrix.rows(); ++r)
    for (size_t c = 0; c < matrix.cols(); ++c)
      out.push_back(matrix[r][c]);
  return out;
}

} // namespace

std::string dump_quantized_llr(const Matrix<float>& llr_mat,
                               const DecodeRequest& request,
                               bool enable,
                               const std::string& path_override,
                               const std::string& default_suffix) {
  if (!enable) return {};

  std::filesystem::path out_path =
      path_override.empty()
          ? std::filesystem::path("data/llr/quantized_llr_" + std::string(request.label) + default_suffix + ".txt")
          : std::filesystem::path(path_override);

  if (out_path.has_parent_path()) {
    std::error_code ec;
    std::filesystem::create_directories(out_path.parent_path(), ec);
    if (ec && !request.quiet) {
      std::cout << "[WARN] (" << request.label
                << ") Failed to create dir for quantized LLR: "
                << ec.message() << "\n";
    }
  }

  std::ofstream ofs(out_path, std::ios::out | std::ios::trunc);
  if (!ofs.is_open()) {
    if (!request.quiet) {
      std::cout << "[WARN] (" << request.label
                << ") Cannot open file to dump quantized LLR: "
                << out_path.string() << "\n";
    }
    return {};
  }

  for (float v : flatten_matrix(llr_mat)) {
    ofs << v << '\n';
  }

  if (!request.quiet) {
    std::cout << "[INFO] (" << request.label << ") Quantized LLR saved to "
              << out_path.string() << " (" << llr_mat.rows() * llr_mat.cols()
              << " values)\n";
    std::cout << "        MATLAB: x = readmatrix('" << out_path.string()
              << "'); histogram(x);\n";
  }

  return out_path.string();
}

} // namespace newcode
