#pragma once

#include <string>
#include <vector>

#include "new_float_only/params.hpp"
#include "new_float_only/rx/ber/ber.hpp"

namespace new_float_only {

/**
 * 单次运行的最终结果。
 * 该结果面向实验脚本或示例程序，包含 BER、错误位置以及日志/矩阵输出路径。
 */
struct SingleRunResult {
  float ebn0_db = 0.0f;
  BerStats pre_fec;
  BerStats post_fec;
  std::vector<std::size_t> pre_fec_error_positions;
  std::vector<std::size_t> post_fec_error_positions;
  std::string log_path;
  std::string work_llr_path;
};

/**
 * 执行一次完整的浮点 oFEC 链路实验。
 * 函数会构建日志文件、运行 pipeline，并把结果回传给调用方。
 */
SingleRunResult run_single(const SingleRunConfig& config);

}  // namespace new_float_only
