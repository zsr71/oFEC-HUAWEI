#pragma once

#include "ofec_tile_input.ipp"
#include "ofec_tile_extrinsic_normalize.ipp"

namespace newcode {
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

inline void log_target_history(const Params& params,
                               const Matrix<float>& lout,
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

template <typename LLR>
DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>
decode_tile(const TilePrepared<LLR>& prep,
            bool use_hard_decode,
            bool normalize_extrinsic,
            const Params& p,
            const std::vector<bool>* early_stop_row_flags,
            CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn)
{
  auto decoder_res = core_fn(prep.lin_matrix,
                             prep.lch_matrix,
                             use_hard_decode,
                             prep.params_for_core,
                             early_stop_row_flags);

  if (normalize_extrinsic && !use_hard_decode)
  {
    normalize_extrinsic_lout(decoder_res.lout, decoder_res.produced_rows, p.beta);
  }
  {
    // Apply global α scaling on extrinsic outputs (moved from Chase decoder).
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
    // Quantize extrinsics back to the target LLR precision before writing out.
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

} // namespace detail
} // namespace newcode
