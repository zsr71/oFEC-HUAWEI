#pragma once

#include "ofec_tile_input.ipp"
#include "ofec_tile_extrinsic_normalize.ipp"

#include <cstdint>

namespace newcode {
namespace detail {

inline std::string sanitize_label(const std::string& label) {
  // 输入:
  // - label: 追踪目标的原始标签。
  // 输出:
  // - 适合作为文件名的一段安全字符串。
  // 用途:
  // - Chase 调试导出 CSV 时，避免标签中出现路径或特殊字符。
  std::string safe = label.empty() ? std::string("target") : label;
  for (char& ch : safe) {
    if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '-') {
      ch = '_';
    }
  }
  return safe;
}

inline std::string quantized_code_csv_cell(float value,
                                           const newcode::Params& params) {
  // 输入:
  // - value: 已经经过公共量化/裁剪后的 extrinsic 浮点值。
  // - params: 当前 tile 参数，提供 LLR 位宽和 clip。
  // 输出:
  // - 若当前是 qfloat 路径，则返回对应的内部 code；否则返回空串。
  // 用途:
  // - 在 target_*.csv 里补充“这一轮实际写回矩阵前，对应的量化码值”。
  const float clip = params.LLR_CLIP;
  switch (params.LLR_BITS) {
    case 2:  return std::to_string(qfloat::qfloat<2>::from_float(value, clip).code());
    case 3:  return std::to_string(qfloat::qfloat<3>::from_float(value, clip).code());
    case 4:  return std::to_string(qfloat::qfloat<4>::from_float(value, clip).code());
    case 5:  return std::to_string(qfloat::qfloat<5>::from_float(value, clip).code());
    case 6:  return std::to_string(qfloat::qfloat<6>::from_float(value, clip).code());
    case 7:  return std::to_string(qfloat::qfloat<7>::from_float(value, clip).code());
    case 8:  return std::to_string(qfloat::qfloat<8>::from_float(value, clip).code());
    case 9:  return std::to_string(qfloat::qfloat<9>::from_float(value, clip).code());
    case 10: return std::to_string(qfloat::qfloat<10>::from_float(value, clip).code());
    case 11: return std::to_string(qfloat::qfloat<11>::from_float(value, clip).code());
    case 12: return std::to_string(qfloat::qfloat<12>::from_float(value, clip).code());
    case 13: return std::to_string(qfloat::qfloat<13>::from_float(value, clip).code());
    case 14: return std::to_string(qfloat::qfloat<14>::from_float(value, clip).code());
    case 15: return std::to_string(qfloat::qfloat<15>::from_float(value, clip).code());
    default: return std::string();
  }
}

inline std::string csv_float_cell(bool has_value, float value) {
  return has_value ? std::to_string(value) : std::string();
}

inline std::string csv_int_cell(bool has_value, int value) {
  return has_value ? std::to_string(value) : std::string();
}

inline void log_target_history(const newcode::Params& params,
                               const matrix::Matrix<float>& lout,
                               const std::vector<bool>& produced_rows) {
  // 输入:
  // - params: 调试配置来源。
  // - lout: 当前 tile core 输出的外信息矩阵。
  // - produced_rows: 每个 decoder row 是否真的产出结果。
  // 输出:
  // - 无返回值；满足条件时向 CSV 文件追加一行追踪记录。
  // 用途:
  // - 将指定坐标/目标比特在每次 Chase 调用中的 extrinsic 变化落盘，便于离线分析。
  const auto& trace = params.debug_trace;
  if (!trace.enable || !trace.dump_chase_csv) return;
  if (trace.active_chase_entries.empty()) return;
  if (trace.chase_tile_index < 0 || trace.chase_invocation < 0) return;
  namespace fs = std::filesystem;
  const fs::path dir = trace.chase_csv_dir.empty()
                           ? fs::path("data/chase_csv")
                           : fs::path(trace.chase_csv_dir);
  std::error_code ec;
  fs::create_directories(dir, ec);
  const std::string csv_header =
      "invocation,tile_index,global_row,global_col,lin_row_index,lin_k,channel_llr_float,channel_llr_code,extrinsic,extrinsic_code,expected_bit";
  for (const auto& entry : trace.active_chase_entries) {
    if (entry.row_index < 0 || entry.k < 0) continue;
    if (entry.row_index >= static_cast<int>(lout.rows())) continue;
    if (static_cast<size_t>(entry.k) >= lout.cols()) continue;
    if (entry.row_index >= static_cast<int>(produced_rows.size())) continue;
    if (!produced_rows[static_cast<size_t>(entry.row_index)]) continue;
    float extrinsic = lout[static_cast<size_t>(entry.row_index)]
                         [static_cast<size_t>(entry.k)];
    std::string label = entry.label;
    if (label.empty()) {
      if (entry.bit_index >= 0) {
        label = "bit" + std::to_string(entry.bit_index);
      } else {
        label = "row" + std::to_string(entry.global_row) + "_col" +
                std::to_string(entry.global_col);
      }
    }
    const std::string safe_label = sanitize_label(label);
    const fs::path file =
        dir / ("target_" + safe_label + ".csv");
    const bool existed = fs::exists(file);
    bool append_mode = existed;
    if (existed) {
      std::ifstream in(file);
      std::string first_line;
      if (!std::getline(in, first_line) || first_line != csv_header) {
        // 旧格式 target.csv 缺少新列时，直接用新表头重建，避免新旧列数混杂。
        append_mode = false;
      }
    }
    std::ofstream out(file, append_mode ? std::ios::app : std::ios::trunc);
    if (!out) continue;
    if (!append_mode) {
      out << csv_header << '\n';
    }
    out << trace.chase_invocation << ',' << trace.chase_tile_index << ','
        << entry.global_row << ',' << entry.global_col << ','
        << entry.row_index << ',' << entry.k << ','
        << csv_float_cell(entry.has_channel_llr_float, entry.channel_llr_float) << ','
        << csv_int_cell(entry.has_channel_llr_code, entry.channel_llr_code) << ','
        << extrinsic << ',' << quantized_code_csv_cell(extrinsic, params)
        << ',' << entry.expected_bit << '\n';
  }
}

template <typename LLR>
chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>
decode_tile(const TilePrepared<LLR>& prep,
            bool use_hard_decode,
            bool normalize_extrinsic,
            const newcode::Params& p,
            const std::vector<bool>* early_stop_row_flags,
            const std::vector<uint8_t>* mux_state,
            CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn)
{
  // 输入:
  // - prep: 已经排布好的 tile 输入，包括 lin/lch 矩阵和 trace 信息。
  // - use_hard_decode: 是否切换到硬判回退。
  // - normalize_extrinsic: 是否对 soft 路径外信息做归一化。
  // - p: 当前 tile 参数。
  // - early_stop_row_flags: 每个 decoder row 的 early-stop 标志。
  // - mux_state: 当前 tile 的 mux 分配状态。
  // - core_fn: 实际执行 Chase/硬判回退的 core。
  // 输出:
  // - DecoderCoreResult，包含 lout 和 produced_rows。
  // 用途:
  // - 这是 tile 级真正调用 decoder core 的地方，并负责对 core 输出做后处理。
  auto decoder_res = core_fn(prep.lin_matrix,
                             prep.lch_matrix,
                             use_hard_decode,
                             prep.params_for_core,
                             early_stop_row_flags,
                             mux_state);

  if (normalize_extrinsic && !use_hard_decode)
  {
    // 仅 soft 路径做归一化；硬判回退输出保持固定幅度语义。
    normalize_extrinsic_lout(decoder_res.lout, decoder_res.produced_rows, p.beta);
  }
  {
    // 对外部输出应用全局 α 缩放（从 Chase 解码器移出）。
    const std::size_t Rcnt = decoder_res.lout.rows();
    const std::size_t Ccnt = decoder_res.lout.cols();
    for (std::size_t r = 0; r < Rcnt; ++r)
    {
      if (!decoder_res.produced_rows[r]) continue;
      for (std::size_t j = 0; j < Ccnt; ++j)
      {
        decoder_res.lout[r][j] *= p.ALPHA;
      }
    }
  }

  {
    // 在写出之前将外信息量化/裁剪回目标 LLR 精度。
    const std::size_t Rcnt = decoder_res.lout.rows();
    const std::size_t Ccnt = decoder_res.lout.cols();
    for (std::size_t r = 0; r < Rcnt; ++r)
    {
      if (!decoder_res.produced_rows[r]) continue;
      for (std::size_t j = 0; j < Ccnt; ++j)
      {
        decoder_res.lout[r][j] =
            ExtrinsicQuantizer<LLR>::quantize(decoder_res.lout[r][j]);
      }
    }
  }

  log_target_history(prep.params_for_core, decoder_res.lout,
                     decoder_res.produced_rows);

  return decoder_res;
}

} // 命名空间 detail
} // 命名空间 newcode
