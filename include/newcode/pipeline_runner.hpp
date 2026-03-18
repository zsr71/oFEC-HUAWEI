#pragma once

#include <limits>
#include <string>
#include <vector>

#include "newcode/rx/ber/ber.hpp"
#include "newcode/decoder_api.hpp"
#include "newcode/params.hpp"

namespace newcode {

struct PipelineConfig {
  std::string decoder_name = "chase_baseline";
  std::string interleaver_name = "ofec";
  bool normalize_extrinsic = true;
  unsigned bits_per_symbol = 2;
  bool quiet = false;
  bool dump_quantized_llr = false;           // 打开时保存量化后的 LLR（便于画直方图）
  std::string quantized_llr_output_path;     // 输出路径（为空则由调用方确定）
  bool dump_work_llr = false;                // 是否保存 decoder 工作矩阵 work_llr
  std::string work_llr_output_path;          // work_llr 输出路径（为空则由调用方确定）
};

struct PipelineResult {
  float ebn0_db = std::numeric_limits<float>::quiet_NaN();
  BerStats pre_fec;
  BerStats post_fec;
  std::vector<std::size_t> pre_fec_error_positions;
  std::vector<std::size_t> post_fec_error_positions;
  std::vector<double> tile_early_stop_pct;
  std::vector<double> tile_row_early_stop_pct;
  std::string dequantized_llr_path; // 反量化后的 LLR 保存位置
  std::string float_llr_path;       // 解调 float LLR 的保存位置
  std::string quantized_codes_path; // 量化码字的保存位置
};

inline constexpr float DEFAULT_EBN0_DB = 3.24f;

PipelineResult run_pipeline(const Params& params,
                            const std::string& label,
                            float ebn0_dB = DEFAULT_EBN0_DB);

PipelineResult run_pipeline(const Params& params,
                            const PipelineConfig& config,
                            const std::string& label,
                            float ebn0_dB = DEFAULT_EBN0_DB);

} // namespace newcode
