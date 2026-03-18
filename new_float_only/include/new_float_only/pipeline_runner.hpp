#pragma once

#include <limits>
#include <string>
#include <vector>

#include "new_float_only/params.hpp"
#include "new_float_only/rx/ber/ber.hpp"

namespace new_float_only {

/**
 * 完整 TX->Channel->RX 单次链路的运行配置。
 * 这里只放链路级参数；解码细节仍由 DecoderConfig 提供。
 */
struct PipelineConfig {
  bool normalize_extrinsic = true;
  unsigned bits_per_symbol = 2;
  int bitgen_seed = 56456;
  int channel_seed = 57112;
  bool generate_random_bits = true;
  bool collect_error_positions = true;
  bool quiet = false;
};

/**
 * 单次链路运行结果。
 * 除了 pre/post BER 之外，还会保留错误位置和可选输出文件路径。
 */
struct PipelineResult {
  float ebn0_db = std::numeric_limits<float>::quiet_NaN();
  BerStats pre_fec;
  BerStats post_fec;
  // 当 collect_error_positions=false 时，这两个数组会保持为空，用于减少 sweep 场景下的内存和比较开销。
  std::vector<std::size_t> pre_fec_error_positions;
  std::vector<std::size_t> post_fec_error_positions;
  std::string work_llr_path;
};

PipelineResult run_pipeline(const Params& params,
                            const PipelineConfig& config,
                            const std::string& label,
                            float ebn0_dB);

}  // namespace new_float_only
