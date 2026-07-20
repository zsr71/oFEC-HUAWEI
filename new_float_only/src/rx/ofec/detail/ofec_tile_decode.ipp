#pragma once

#include "ofec_tile_input.ipp"
#include "ofec_tile_extrinsic_normalize.ipp"

#include <cstdint>

namespace new_float_only {
namespace detail {

inline std::string sanitize_label(const std::string& label) {
  std::string safe = label.empty() ? std::string("target") : label;
  for (char& ch : safe) {
    if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '-') {
      ch = '_';
    }
  }
  return safe;
}

inline void log_target_history(const new_float_only::Params& params,
                               const matrix::Matrix<float>& lout,
                               const std::vector<bool>& produced_rows) {
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
    std::ofstream out(file, std::ios::app);
    if (!out) continue;
    if (!existed) {
      out << "invocation,tile_index,global_row,global_col,lin_row_index,lin_k,extrinsic,expected_bit\n";
    }
    out << trace.chase_invocation << ',' << trace.chase_tile_index << ','
        << entry.global_row << ',' << entry.global_col << ','
        << entry.row_index << ',' << entry.k << ','
        << extrinsic << ',' << entry.expected_bit << '\n';
  }
}

/**
 * 对当前 Tile 的重排结果执行行级 core，并完成 extrinsic 后处理。
 * 关键语句：
 * 1. normalize_extrinsic_lout 用于统一外信息幅度；
 * 2. lout *= ALPHA 对应论文里的外信息缩放；
 * 3. 这里只返回 extrinsic，不在此处叠加信道项。
 */
chase::DecoderCoreResult
decode_tile(const TilePrepared& prep,
            bool use_hard_decode,
            bool normalize_extrinsic,
            const new_float_only::Params& p,
            CoreFn core_fn)
{
  auto decoder_res = core_fn(prep.lin_matrix,
                             prep.lch_matrix,
                             use_hard_decode,
                             prep.params_for_core);

  if (normalize_extrinsic && !use_hard_decode)
  {
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
    // 这里不做量化回写，新库固定使用 float 通路。
    const std::size_t Rcnt = decoder_res.lout.rows();
    const std::size_t Ccnt = decoder_res.lout.cols();
    for (std::size_t r = 0; r < Rcnt; ++r)
    {
      if (!decoder_res.produced_rows[r]) continue;
      for (std::size_t j = 0; j < Ccnt; ++j)
      {
        decoder_res.lout[r][j] = decoder_res.lout[r][j];
      }
    }
  }

  log_target_history(prep.params_for_core, decoder_res.lout,
                     decoder_res.produced_rows);

  return decoder_res;
}

} // 命名空间 detail
} // 命名空间 newcode
