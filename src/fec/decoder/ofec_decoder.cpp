#include "newcode/ofec_decoder.hpp"
#include "newcode/chase256.hpp" // 保留 Chase 头；本文档内有三参前向声明
#include "newcode/ofec_decoder_hard.hpp"
#include "newcode/decoder_core.hpp"
#include "newcode/llr_utils.hpp"
#include "newcode/decoder_api.hpp"
#include "newcode/quantized_llr_dump.hpp"
#include "common/lin_matrix_utils.hpp"

#include <filesystem>
#include <fstream>
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <cctype>
#include <memory>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace newcode {

namespace {

template <typename LLR>
using CoreFn = DecoderCoreResult<LLR> (*)(const Matrix<LLR>&,
                                          const Matrix<LLR>&,
                                          bool,
                                          const Params&);

template <typename LLR>
TileProcessResult<LLR> process_tile_impl(const Matrix<LLR>& tile_in,
                                         const Matrix<LLR>& ch_tile,
                                         const Params& p,
                                         size_t tile_top_row_global,
                                         bool use_hard_decode,
                                         bool normalize_extrinsic,
                                         const Matrix<float>* tx_llr_ref,
                                         CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
                                         Matrix<float>* last_tile_history_accum,
                                         bool capture_last_tile_history)
{
  using Adapter  = LinMatrixAdapter<LLR>;
  using CoreLLR  = typename Adapter::core_type;
  auto core_to_float = [](const CoreLLR& value) -> float {
    if constexpr (std::is_same_v<CoreLLR, float> || std::is_same_v<CoreLLR, double>) {
      return static_cast<float>(value);
    } else {
      return llr_to_float(value);
    }
  };

  constexpr int B         = static_cast<int>(Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(Params::NUM_SUBBLOCK_COLS * B);            // 128
  constexpr int K         = static_cast<int>(Params::BCH_K);                             // 239
  constexpr int TAKE_BITS = K - N;                                                       // 111
  constexpr int BCH_PAR   = static_cast<int>(Params::BCH_PARITY_BITS);                   // 16
  constexpr int OVR_IDX   = static_cast<int>(Params::BCH_OVERALL_IDX);                   // 255

  const size_t H = tile_in.rows();
  const size_t W = tile_in.cols();
  assert(W == static_cast<size_t>(N));
  assert(ch_tile.rows() == H && ch_tile.cols() == W);

  Matrix<LLR> tile_out = tile_in;

  const auto& trace_cfg = p.debug_trace;
  const bool trace_has_coords = (trace_cfg.row >= 0 && trace_cfg.col >= 0);
  const bool trace_has_targets = !trace_cfg.targets.empty();
  const bool trace_enabled = trace_cfg.enable && (trace_has_coords || trace_has_targets);
  const long trace_row = trace_has_coords ? trace_cfg.row : -1;
  const long trace_col = trace_has_coords ? trace_cfg.col : -1;
  Params params_for_core = p;  // may be updated
  params_for_core.debug_trace.active_chase_entries.clear();
  if (params_for_core.debug_trace.chase_expected_bits) {
    params_for_core.debug_trace.chase_expected_bits.reset();
    params_for_core.debug_trace.chase_expected_bits_row = nullptr;
  }
  std::shared_ptr<std::vector<std::vector<int8_t>>> expected_bits;
  auto should_trace = [&](bool flag, long rr, long cc) -> bool {
    return trace_enabled && flag && rr == trace_row && cc == trace_col;
  };
  auto register_chase_entry = [&](const std::string& label,
                                  long bit_index,
                                  long rr, long cc,
                                  size_t row_idx, int k) {
    int expected_bit = -1;
    if (expected_bits && row_idx < expected_bits->size() &&
        k >= 0 && static_cast<size_t>(k) < (*expected_bits)[row_idx].size()) {
      expected_bit = (*expected_bits)[row_idx][static_cast<size_t>(k)];
    }
    auto& entries = params_for_core.debug_trace.active_chase_entries;
    for (const auto& entry : entries) {
      if (entry.row_index == static_cast<int>(row_idx) &&
          entry.k == k &&
          entry.global_row == rr &&
          entry.global_col == cc &&
          entry.bit_index == bit_index &&
          entry.label == label &&
          entry.expected_bit == expected_bit) {
        return;
      }
    }
    newcode::Params::DebugTraceConfig::ChaseTraceEntry entry;
    entry.row_index = static_cast<int>(row_idx);
    entry.k = k;
    entry.global_row = rr;
    entry.global_col = cc;
    entry.bit_index = bit_index;
    entry.label = label;
    entry.expected_bit = expected_bit;
    params_for_core.debug_trace.active_chase_entries.push_back(entry);
  };
  auto try_mark_chase_coord = [&](long rr, long cc, size_t row_idx, int k) {
    if (!trace_enabled) return;
    if (trace_has_coords && rr == trace_row && cc == trace_col) {
      const std::string label =
          "row" + std::to_string(rr) + "_col" + std::to_string(cc);
      register_chase_entry(label, -1, rr, cc, row_idx, k);
    }
    for (size_t idx = 0; idx < trace_cfg.targets.size(); ++idx) {
      const auto& target = trace_cfg.targets[idx];
      if (target.row == rr && target.col == cc) {
        std::string label;
        if (!target.label.empty()) {
          label = target.label;
        } else if (target.bit_index >= 0) {
          label = "bit" + std::to_string(target.bit_index);
        } else {
          label = "target" + std::to_string(idx);
        }
        register_chase_entry(label, target.bit_index, rr, cc, row_idx, k);
      }
    }
  };
  auto sanitize_label = [](const std::string& label) {
    std::string safe = label.empty() ? std::string("target") : label;
    for (char& ch : safe) {
      if (!std::isalnum(static_cast<unsigned char>(ch)) && ch != '-') {
        ch = '_';
      }
    }
    return safe;
  };
  auto log_target_history = [&](const Params& params,
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
  };

  const int SBR = p.CHASE_SBR;
  if (SBR != 1 && SBR != 2)
      throw std::invalid_argument("process_tile: CHASE_SBR must be 1 or 2.");

  const size_t rows_to_decode = static_cast<size_t>(SBR) * static_cast<size_t>(B);
  if (rows_to_decode == 0) {
      return TileProcessResult<LLR>{tile_out, false};
  }

  const size_t decoder_cols = static_cast<size_t>(2 * N); // 256
  if (tx_llr_ref) {
    expected_bits = std::make_shared<std::vector<std::vector<int8_t>>>(
        rows_to_decode, std::vector<int8_t>(decoder_cols, -1));
  }
  Matrix<CoreLLR> lin_matrix(rows_to_decode, decoder_cols);
  Matrix<CoreLLR> lch_matrix(rows_to_decode, decoder_cols);
  std::vector<size_t> row_local_lookup(rows_to_decode, 0);
  std::vector<size_t> row_global_lookup(rows_to_decode, 0);

  for (int s = 0; s < SBR; ++s){
      const size_t sbr_row0_local = H - static_cast<size_t>((SBR-s) * B);
      for (int r_off = 0; r_off < B; ++r_off)
      {
          const size_t row_idx   = static_cast<size_t>(s * B + r_off);
          const size_t row_local = sbr_row0_local + static_cast<size_t>(r_off);
          const size_t row_global = tile_top_row_global + row_local;

          row_local_lookup[row_idx]  = row_local;
          row_global_lookup[row_idx] = row_global;

          const long R = static_cast<long>(row_global / static_cast<size_t>(B));
          const int  r = static_cast<int>(row_global % static_cast<size_t>(B));

          for (int k = 0; k < N; ++k)
          {
              const long br = (R ^ 1L)
                            - static_cast<long>(2 * p.NUM_GUARD_SUBROWS)
                            - static_cast<long>(2 * (N / B))
                            + static_cast<long>(2 * (k / B));
              const long bc = static_cast<long>(k / B);
              const long bit_row_in_block = static_cast<long>((k % B) ^ r);
              const long bit_col_in_block = static_cast<long>(r);

              const long rr_global = br * B + bit_row_in_block;
              const long cc_global = bc * B + bit_col_in_block;

              if (should_trace(trace_cfg.log_read_mapping, rr_global, cc_global)) {
                  std::cout << " AS History READ Mapping k=" << k
                            << " to global pos (" << rr_global << "," << cc_global << ")" << '\n';
              }
              try_mark_chase_coord(rr_global, cc_global, row_idx, k);



              const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
              const long cc_local2 = cc_global;

              if (rr_local2 >= 0 && rr_local2 < static_cast<long>(H) &&
                  cc_local2 >= 0 && cc_local2 < static_cast<long>(W))
              {
              const LLR Lch = ch_tile[static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)];
              const LLR La  = tile_in [static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)];
              lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
              lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
              }
              else {
                  throw std::out_of_range("process_tile: old info position out of tile range.");
              }
              if (expected_bits &&
                  rr_global >= 0 && cc_global >= 0 &&
                  static_cast<size_t>(rr_global) < tx_llr_ref->rows() &&
                  static_cast<size_t>(cc_global) < tx_llr_ref->cols()) {
                  const float v = (*tx_llr_ref)[static_cast<size_t>(rr_global)]
                                              [static_cast<size_t>(cc_global)];
                  (*expected_bits)[row_idx][static_cast<size_t>(k)] = (v < 0.0f) ? 1 : 0;
              }
          }

          for (int i = 0; i < TAKE_BITS; ++i) {
              const int k = N + i;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t src_col = Ct * static_cast<size_t>(B) + ct;
              if (should_trace(trace_cfg.log_read_mapping, row_local+tile_top_row_global, src_col)) {
                  std::cout << " AS New Information READ Mapping k=" << k
                            << " to global pos (" << (row_local + tile_top_row_global) << "," << src_col << ")" << '\n';
              }
              try_mark_chase_coord(static_cast<long>(row_local + tile_top_row_global),
                                   static_cast<long>(src_col), row_idx, k);
              const LLR Lch = ch_tile[row_local][src_col];
              const LLR La  = tile_in [row_local][src_col];
              lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
              lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
              if (expected_bits &&
                  (tile_top_row_global + row_local) < tx_llr_ref->rows() &&
                  src_col < tx_llr_ref->cols()) {
                  const float v = (*tx_llr_ref)[tile_top_row_global + row_local][src_col];
                  (*expected_bits)[row_idx][static_cast<size_t>(k)] = (v < 0.0f) ? 1 : 0;
              }
          }

          for (int j = 0; j < BCH_PAR; ++j) {
              const int k = K + j;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t src_col = Ct * static_cast<size_t>(B) + ct;

              const LLR Lch = ch_tile[row_local][src_col];
              const LLR La  = tile_in [row_local][src_col];
              lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
              lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
              if (expected_bits &&
                  (tile_top_row_global + row_local) < tx_llr_ref->rows() &&
                  src_col < tx_llr_ref->cols()) {
                  const float v = (*tx_llr_ref)[tile_top_row_global + row_local][src_col];
                  (*expected_bits)[row_idx][static_cast<size_t>(k)] = (v < 0.0f) ? 1 : 0;
              }
          }

          {
              const int k = OVR_IDX;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t src_col = Ct * static_cast<size_t>(B) + ct;

              const LLR Lch = ch_tile[row_local][src_col];
              const LLR La  = tile_in [row_local][src_col];
              lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
              lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
              if (expected_bits &&
                  (tile_top_row_global + row_local) < tx_llr_ref->rows() &&
                  src_col < tx_llr_ref->cols()) {
                  const float v = (*tx_llr_ref)[tile_top_row_global + row_local][src_col];
                  (*expected_bits)[row_idx][static_cast<size_t>(k)] = (v < 0.0f) ? 1 : 0;
              }
          }
      }
    }

  if (expected_bits) {
    params_for_core.debug_trace.chase_expected_bits = expected_bits;
  } else {
    params_for_core.debug_trace.chase_expected_bits.reset();
  }
  if (!params_for_core.debug_trace.active_chase_entries.empty()) {
    params_for_core.debug_trace.chase_decoder_row =
        params_for_core.debug_trace.active_chase_entries.front().row_index;
    params_for_core.debug_trace.chase_decoder_col =
        params_for_core.debug_trace.active_chase_entries.front().k;
  } else {
    params_for_core.debug_trace.chase_decoder_row = -1;
    params_for_core.debug_trace.chase_decoder_col = -1;
  }

  bool early_stop_triggered = tile_should_early_stop(lin_matrix);

  auto decoder_res = core_fn(lin_matrix, lch_matrix, use_hard_decode, params_for_core);

  // Normalization only applies to soft extrinsics; skip for hard-decoded tiles.
  if (normalize_extrinsic && !use_hard_decode)
  {
      auto is_fallback = [&](float w) -> bool {
          const float target = p.beta;
          const float diff   = std::fabs(std::fabs(w) - target);
          const float tol    = 1e-4f * std::max(1.0f, target);
          return diff <= tol;
      };

      double acc = 0.0;
      std::size_t cnt = 0;

      const std::size_t Rcnt = decoder_res.lout.rows();
      const std::size_t Ccnt = decoder_res.lout.cols();

      for (std::size_t r = 0; r < Rcnt; ++r)
      {
          if (!decoder_res.produced_rows[r]) continue;
          for (std::size_t j = 0; j < Ccnt; ++j)
          {
              const float w = decoder_res.lout[r][j];
              if (is_fallback(w)) {acc += std::fabs(w/p.beta); ++cnt; continue;};
              acc += std::fabs(w);
              ++cnt;
          }
      }

      if (cnt > 0)
      {
          const float g_alpha = static_cast<float>(acc / static_cast<double>(cnt));
          if (g_alpha > 0.f)
          {
              const float scale = 1.0f / g_alpha;
              for (std::size_t r = 0; r < Rcnt; ++r)
              {
                  if (!decoder_res.produced_rows[r]) continue;
                  for (std::size_t j = 0; j < Ccnt; ++j)
                  {
                      //if (is_fallback(decoder_res.lout[r][j])) continue;
                      decoder_res.lout[r][j] *= scale;
                  }
              }
          }
      }
  }
  if(1)
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

  log_target_history(params_for_core, decoder_res.lout,
                     decoder_res.produced_rows);

  for (int s = 0; s < SBR; ++s)
  {
      for (int r_off = 0; r_off < B; ++r_off)
      {
          const size_t row_idx = static_cast<size_t>(s * B + r_off);
          if (row_idx >= decoder_res.produced_rows.size()) continue;
          if (!decoder_res.produced_rows[row_idx]) continue;

          const size_t row_local  = row_local_lookup[row_idx];
          const size_t row_global = row_global_lookup[row_idx];
          const auto&  lout_row   = decoder_res.lout[row_idx];

          const int r = static_cast<int>(row_global % static_cast<size_t>(B));

          for (int i = 0; i < TAKE_BITS; ++i) {
              const int k = N + i;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t col = Ct * static_cast<size_t>(B) + ct;
              const LLR prior_llr = tile_out[row_local][col];
              const LLR extrinsic_llr =
                  llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
              tile_out[row_local][col] = extrinsic_llr;

              if (capture_last_tile_history && last_tile_history_accum) {
                  const long rr_global = static_cast<long>(row_global);
                  const long cc_global = static_cast<long>(col);
                  if (rr_global >= 0 && cc_global >= 0) {
                      const size_t rr_idx_global = static_cast<size_t>(rr_global);
                      const size_t cc_idx_global = static_cast<size_t>(cc_global);
                      if (rr_idx_global < last_tile_history_accum->rows() &&
                          cc_idx_global < last_tile_history_accum->cols()) {
                          const CoreLLR combined =
                              Adapter::combine(extrinsic_llr, prior_llr);
                          (*last_tile_history_accum)[rr_idx_global][cc_idx_global] =
                              core_to_float(combined);
                      }
                  }
              }
          }
          for (int j = 0; j < BCH_PAR; ++j) {
              const int k = K + j;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t col = Ct * static_cast<size_t>(B) + ct;
              const LLR prior_llr = tile_out[row_local][col];
              const LLR extrinsic_llr =
                  llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
              tile_out[row_local][col] = extrinsic_llr;

              if (capture_last_tile_history && last_tile_history_accum) {
                  const long rr_global = static_cast<long>(row_global);
                  const long cc_global = static_cast<long>(col);
                  if (rr_global >= 0 && cc_global >= 0) {
                      const size_t rr_idx_global = static_cast<size_t>(rr_global);
                      const size_t cc_idx_global = static_cast<size_t>(cc_global);
                      if (rr_idx_global < last_tile_history_accum->rows() &&
                          cc_idx_global < last_tile_history_accum->cols()) {
                          const CoreLLR combined =
                              Adapter::combine(extrinsic_llr, prior_llr);
                          (*last_tile_history_accum)[rr_idx_global][cc_idx_global] =
                              core_to_float(combined);
                      }
                  }
              }
          }
          {
              const int k = OVR_IDX;
              const size_t Ct = static_cast<size_t>((k - N) / B);
              const size_t ct = static_cast<size_t>((k % B) ^ r);
              const size_t col = Ct * static_cast<size_t>(B) + ct;
              const LLR prior_llr = tile_out[row_local][col];
              const LLR extrinsic_llr =
                  llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
              tile_out[row_local][col] = extrinsic_llr;

              if (capture_last_tile_history && last_tile_history_accum) {
                  const long rr_global = static_cast<long>(row_global);
                  const long cc_global = static_cast<long>(col);
                  if (rr_global >= 0 && cc_global >= 0) {
                      const size_t rr_idx_global = static_cast<size_t>(rr_global);
                      const size_t cc_idx_global = static_cast<size_t>(cc_global);
                      if (rr_idx_global < last_tile_history_accum->rows() &&
                          cc_idx_global < last_tile_history_accum->cols()) {
                          const CoreLLR combined =
                              Adapter::combine(extrinsic_llr, prior_llr);
                          (*last_tile_history_accum)[rr_idx_global][cc_idx_global] =
                              core_to_float(combined);
                      }
                  }
              }
          }

          const long R = static_cast<long>(row_global / static_cast<size_t>(B));
          for (int k = 0; k < N; ++k)
          {
              const long br = (R ^ 1L)
                            - static_cast<long>(2 * p.NUM_GUARD_SUBROWS)
                            - static_cast<long>(2 * (N / B))
                            + static_cast<long>(2 * (k / B));
              const long bc = static_cast<long>(k / B);
              const long bit_row_in_block = static_cast<long>((k % B) ^ r);
              const long bit_col_in_block = static_cast<long>(r);

              const long rr_global = br * B + bit_row_in_block;
              const long cc_global = bc * B + bit_col_in_block;

              if (should_trace(trace_cfg.log_write_mapping, rr_global, cc_global)) {
                  std::cout << " WRITE Mapping k=" << k
                            << " to global pos (" << rr_global << "," << cc_global << ")" << '\n';
                  std::cout << "  Value=" << llr_to_float(lout_row[static_cast<size_t>(k)]) << '\n';
              }


              const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
              const long cc_local2 = cc_global;

              const bool in_range =
                  (rr_local2 >= 0 && rr_local2 < static_cast<long>(H) &&
                   cc_local2 >= 0 && cc_local2 < static_cast<long>(W));

              assert(in_range && "process_tile: write-back out of tile range");

              const size_t rr_idx_local = static_cast<size_t>(rr_local2);
              const size_t cc_idx_local = static_cast<size_t>(cc_local2);
              const LLR prior_llr =
                  tile_out[rr_idx_local][cc_idx_local];
              const LLR extrinsic_llr =
                  llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
              tile_out[rr_idx_local][cc_idx_local] = extrinsic_llr;

              if (capture_last_tile_history && last_tile_history_accum) {
                  if (rr_global >= 0 && cc_global >= 0) {
                      const size_t rr_idx_global = static_cast<size_t>(rr_global);
                      const size_t cc_idx_global = static_cast<size_t>(cc_global);
                      if (rr_idx_global < last_tile_history_accum->rows() &&
                          cc_idx_global < last_tile_history_accum->cols()) {
                          const CoreLLR combined = Adapter::combine(extrinsic_llr, prior_llr);
                          (*last_tile_history_accum)[rr_idx_global][cc_idx_global] =
                              core_to_float(combined);
                      }
                  }
              }
          }
      }
  }

  return TileProcessResult<LLR>{std::move(tile_out), early_stop_triggered};
}

template <typename LLR>
void process_window_impl(Matrix<LLR>& work_llr,
                         const Matrix<LLR>& channel_llr,
                         size_t win_start, size_t win_end, const Params& p,
                         size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                         std::vector<TileEarlyStopCounter>* tile_stats,
                         bool normalize_extrinsic,
                         const Matrix<float>* tx_llr_ref,
                         CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
                         Matrix<float>* last_tile_history_accum)
{
  (void)win_start;
  const auto& trace_cfg = p.debug_trace;
  const bool trace_has_coords = (trace_cfg.row >= 0 && trace_cfg.col >= 0);
  const bool trace_mismatch =
      trace_cfg.enable && trace_cfg.log_mismatch && trace_has_coords;
  const long trace_row = trace_has_coords ? trace_cfg.row : -1;
  const long trace_col = trace_has_coords ? trace_cfg.col : -1;
  const size_t N = Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;
  static size_t chase_invocation_counter = 0;

  auto pick_float = [](const std::vector<float>& tbl, size_t idx, float fallback) -> float {
      return (idx < tbl.size()) ? tbl[idx] : fallback;
  };
  auto pick_int = [](const std::vector<int>& tbl, size_t idx, int fallback) -> int {
      return (idx < tbl.size()) ? tbl[idx] : fallback;
  };
  std::vector<bool> hard_tile_mask(TILES_PER_WIN);
  int last_soft_tile_idx = -1;
  for (size_t t = 0; t < TILES_PER_WIN; ++t) {
      const bool is_hard = pick_int(p.HARD_TILE_LIST, t, p.HARD_DECODE_DEFAULT ? 1 : 0) != 0;
      hard_tile_mask[t] = is_hard;
      if (!is_hard) {
          last_soft_tile_idx = static_cast<int>(t);
      }
  }

  for (size_t t = 0; t < TILES_PER_WIN; ++t)
  {
        const size_t tile_bottom_row = win_end  - t * tile_stride_rows;
        const size_t tile_top_row    = tile_bottom_row + 1 - tile_height_rows;

        const size_t tile_height_rows_actual = tile_bottom_row - tile_top_row + 1;
        Matrix<LLR> tile_in(tile_height_rows_actual, N);
        Matrix<LLR> ch_tile(tile_height_rows_actual, N);

        const bool use_hard = hard_tile_mask[t];
        const bool use_history_input =
            use_hard && last_tile_history_accum &&
            last_soft_tile_idx >= 0 &&
            static_cast<int>(t) > last_soft_tile_idx;

        for (size_t r = 0; r < tile_height_rows_actual; ++r) {
            for (size_t c = 0; c < N; ++c) {
                const size_t global_row = tile_top_row + r;
                if (use_history_input) {
                    const float hist = (*last_tile_history_accum)[global_row][c];
                    tile_in[r][c] = llr_from_float<LLR>(hist);
                } else {
                    tile_in[r][c]  = work_llr[global_row][c];
                }
                ch_tile[r][c]  = channel_llr[global_row][c];
            }
        }

        Params tile_params = p;
        tile_params.beta = pick_float(p.beta_list, t, p.beta);
        tile_params.ALPHA = pick_float(p.ALPHA_LIST, t, p.ALPHA);
        tile_params.debug_trace.chase_tile_index = static_cast<int>(t);
        tile_params.debug_trace.chase_invocation =
            static_cast<int>(++chase_invocation_counter);

        const bool capture_history =
            last_tile_history_accum &&
            (use_hard || static_cast<int>(t) == last_soft_tile_idx);

        TileProcessResult<LLR> tile_result = process_tile_impl<LLR>(tile_in, ch_tile, tile_params,
                                                                /*tile_top_row_global=*/tile_top_row,
                                                                /*use_hard_decode=*/use_hard,
                                                                /*normalize_extrinsic=*/normalize_extrinsic,
                                                                tx_llr_ref,
                                                                core_fn,
                                                                last_tile_history_accum,
                                                                capture_history);

    if (tile_stats && t < tile_stats->size()) {
      auto& counter = (*tile_stats)[t];
      counter.total += 1;
      if (tile_result.early_stop_triggered) {
        counter.triggered += 1;
      }
    }

    for (size_t r = 0; r < tile_height_rows_actual; ++r) {
      const size_t global_row = tile_top_row + r;
      for (size_t c = 0; c < work_llr.cols(); ++c) {
        const auto incoming = tile_result.tile_out[r][c];
        if (trace_mismatch &&
            static_cast<long>(global_row) == trace_row &&
            static_cast<long>(c) == trace_col) {
          const float existing_val = llr_to_float(work_llr[global_row][c]);
          const float incoming_val = llr_to_float(incoming);
          const float channel_val  = llr_to_float(channel_llr[global_row][c]);
          if (existing_val != incoming_val) {
            std::cout << "Mismatch at work_llr[" << trace_row << "][" << trace_col
                      << "]: tile index " << t
                      << " incoming=" << incoming_val << '\n'
                      << " channel =" << channel_val << '\n';
          }
        }
        work_llr[global_row][c] = incoming;
      }
    }
  }
}

template <typename LLR>
Matrix<LLR> ofec_decode_llr_impl(const Matrix<LLR>& llr_mat, const Params& p,
                                 std::vector<TileEarlyStopCounter>* tile_stats,
                                 bool normalize_extrinsic,
                                 const Matrix<float>* tx_llr_ref,
                                 CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn)
{
  const size_t N = Params::NUM_SUBBLOCK_COLS * Params::BITS_PER_SUBBLOCK_DIM;

  const size_t RROWS = llr_mat.rows();
  const size_t CCOLS = llr_mat.cols();
  if (CCOLS != N)
      throw std::invalid_argument("ofec_decode_llr: llr_mat cols != N.");

  assert(p.valid());

  const size_t TILE_HEIGHT_ROWS = p.tile_height_rows();
  const size_t TILE_STRIDE_ROWS = p.tile_stride_rows();
  const size_t WIN_HEIGHT_ROWS  = p.win_height_rows();
  const size_t POP_PUSH_ROWS    = p.pop_push_rows();
  const size_t TILES_PER_WIN    = p.TILES_PER_WIN;

  
  Matrix<LLR> channel_llr = llr_mat;
  Matrix<LLR> work_llr(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
          work_llr[r][c] = llr_from_float<LLR>(0.0f);
  Matrix<float> last_tile_history_llr(RROWS, N);

  if (RROWS < WIN_HEIGHT_ROWS) {
    if (tile_stats) {
      tile_stats->assign(p.TILES_PER_WIN, TileEarlyStopCounter{});
    }
    return channel_llr;
  }

  size_t win_start     = p.initial_win_start_rows();
  const size_t last_ws = RROWS - WIN_HEIGHT_ROWS;

  std::vector<TileEarlyStopCounter> local_tile_stats;
  std::vector<TileEarlyStopCounter>* stats_ptr = nullptr;
  if (tile_stats) {
    local_tile_stats.assign(p.TILES_PER_WIN, TileEarlyStopCounter{});
    stats_ptr = &local_tile_stats;
  }

  while (win_start <= last_ws) {
    const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;

    process_window_impl<LLR>(work_llr, channel_llr,
                             win_start, win_end, p,
                             TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN,
                             stats_ptr,
                             normalize_extrinsic,
                             tx_llr_ref,
                             core_fn,
                             &last_tile_history_llr);

    win_start += POP_PUSH_ROWS;
  }

  if (tile_stats) {
    *tile_stats = std::move(local_tile_stats);
  }

  Matrix<LLR> out(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r) {
    for (size_t c = 0; c < N; ++c) {
      const float sum = llr_to_float(channel_llr[r][c]) + last_tile_history_llr[r][c];
      out[r][c] = llr_from_float<LLR>(sum);
    }
  }

  // 可选：保存窗口累积后的 work_llr（解码前）
  if (p.DUMP_WORK_LLR) {
    Matrix<float> work_float(RROWS, N);
    for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
        work_float[r][c] = llr_to_float(work_llr[r][c]);

    DecodeRequest dump_req{
        .label = "work_llr",
        .channel_llr = work_float,
        .tx_llr_ref = nullptr,
        .params = p,
        .format = LlrFormat::Float,
        .quant_bits = 16,
        .quant_clip = 0.0f,
        .normalize_extrinsic = true,
        .quiet = false,
        .dump_quantized_llr = true,
        .quantized_llr_output_path = p.WORK_LLR_OUTPUT_PATH.empty()
            ? std::string("data/llr/quantized_llr_work.txt")
            : p.WORK_LLR_OUTPUT_PATH,
        .dump_float_llr = false,
        .float_llr_output_path = {},
        .dump_quantized_codes = false,
        .quantized_codes_output_path = {},
        .dump_work_llr = true,
        .work_llr_output_path = {}
    };
    (void)dump_quantized_llr(work_float, dump_req, true,
                             dump_req.quantized_llr_output_path, "_work");
  }

  return out;
}

} // namespace

template <typename LLR>
TileProcessResult<LLR> process_tile_plain(const Matrix<LLR>& tile_in,
                                          const Matrix<LLR>& ch_tile,
                                          const Params& p,
                                          size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic,
                                          const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                           use_hard_decode, normalize_extrinsic,
                           tx_llr_ref,
                           &Decoder_Core_plain<CoreLLR>,
                           /*last_tile_history_accum=*/nullptr,
                           /*capture_last_tile_history=*/false);
}

template <typename LLR>
TileProcessResult<LLR> process_tile_ebchPF(const Matrix<LLR>& tile_in,
                                           const Matrix<LLR>& ch_tile,
                                           const Params& p,
                                          size_t tile_top_row_global,
                                          bool use_hard_decode,
                                          bool normalize_extrinsic,
                                          const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return process_tile_impl(tile_in, ch_tile, p, tile_top_row_global,
                           use_hard_decode, normalize_extrinsic,
                           tx_llr_ref,
                           &Decoder_Core_ebchPF<CoreLLR>,
                           /*last_tile_history_accum=*/nullptr,
                           /*capture_last_tile_history=*/false);
}

template <typename LLR>
void process_window_plain(Matrix<LLR>& work_llr,
                          const Matrix<LLR>& channel_llr,
                          size_t win_start, size_t win_end, const Params& p,
                          size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                          std::vector<TileEarlyStopCounter>* tile_stats,
                          bool normalize_extrinsic,
                          const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  process_window_impl(work_llr, channel_llr, win_start, win_end, p,
                      tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                      tile_stats, normalize_extrinsic, tx_llr_ref,
                      &Decoder_Core_plain<CoreLLR>,
                      /*last_tile_history_accum=*/nullptr);
}

template <typename LLR>
void process_window_ebchPF(Matrix<LLR>& work_llr,
                           const Matrix<LLR>& channel_llr,
                           size_t win_start, size_t win_end, const Params& p,
                          size_t tile_height_rows, size_t tile_stride_rows, size_t TILES_PER_WIN,
                          std::vector<TileEarlyStopCounter>* tile_stats,
                          bool normalize_extrinsic,
                          const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  process_window_impl(work_llr, channel_llr, win_start, win_end, p,
                      tile_height_rows, tile_stride_rows, TILES_PER_WIN,
                      tile_stats, normalize_extrinsic, tx_llr_ref,
                      &Decoder_Core_ebchPF<CoreLLR>,
                      /*last_tile_history_accum=*/nullptr);
}

template <typename LLR>
Matrix<LLR> ofec_decode_llr_plain(const Matrix<LLR>& llr_mat, const Params& p,
                                  std::vector<TileEarlyStopCounter>* tile_stats,
                                  bool normalize_extrinsic,
                                  const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return ofec_decode_llr_impl(llr_mat, p, tile_stats, normalize_extrinsic,
                              tx_llr_ref,
                              &Decoder_Core_plain<CoreLLR>);
}

template <typename LLR>
Matrix<LLR> ofec_decode_llr_ebchPF(const Matrix<LLR>& llr_mat, const Params& p,
                                   std::vector<TileEarlyStopCounter>* tile_stats,
                                   bool normalize_extrinsic,
                                   const Matrix<float>* tx_llr_ref)
{
  using CoreLLR = typename LinMatrixAdapter<LLR>::core_type;
  return ofec_decode_llr_impl(llr_mat, p, tile_stats, normalize_extrinsic,
                              tx_llr_ref,
                              &Decoder_Core_ebchPF<CoreLLR>);
}

// ===== 显式实例化 =====
template TileProcessResult<float>  process_tile_plain<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool, bool, const Matrix<float>*);
template TileProcessResult<int8_t> process_tile_plain<int8_t>(const Matrix<int8_t>&, const Matrix<int8_t>&, const Params&, size_t, bool, bool, const Matrix<float>*);
template TileProcessResult<float>  process_tile_ebchPF<float >(const Matrix<float>&,  const Matrix<float>&,  const Params&, size_t, bool, bool, const Matrix<float>*);
template TileProcessResult<int8_t> process_tile_ebchPF<int8_t>(const Matrix<int8_t>&, const Matrix<int8_t>&, const Params&, size_t, bool, bool, const Matrix<float>*);

template void process_window_plain<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool,
                                           const Matrix<float>*);
template void process_window_plain<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&, size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool,
                                           const Matrix<float>*);
template void process_window_ebchPF<float >(Matrix<float>&,  const Matrix<float>&,  size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool,
                                           const Matrix<float>*);
template void process_window_ebchPF<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&, size_t, size_t, const Params&,
                                           size_t, size_t, size_t, std::vector<TileEarlyStopCounter>*, bool,
                                           const Matrix<float>*);

template Matrix<float>  ofec_decode_llr_plain<float >(const Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);
template Matrix<int8_t> ofec_decode_llr_plain<int8_t>(const Matrix<int8_t>&, const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);
template Matrix<float>  ofec_decode_llr_ebchPF<float >(const Matrix<float>&,  const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);
template Matrix<int8_t> ofec_decode_llr_ebchPF<int8_t>(const Matrix<int8_t>&, const Params&, std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*);

#define INSTANTIATE_QFLOAT(N) \
template TileProcessResult<newcode::qfloat<N>> process_tile_plain<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    const Params&, size_t, bool, bool, const Matrix<float>*); \
template TileProcessResult<newcode::qfloat<N>> process_tile_ebchPF<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    const Params&, size_t, bool, bool, const Matrix<float>*); \
template void process_window_plain<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    size_t, size_t, const Params&, \
    size_t, size_t, size_t, \
    std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*); \
template void process_window_ebchPF<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    size_t, size_t, const Params&, \
    size_t, size_t, size_t, \
    std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*); \
template Matrix<newcode::qfloat<N>> ofec_decode_llr_plain<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Params&, \
    std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*); \
template Matrix<newcode::qfloat<N>> ofec_decode_llr_ebchPF<newcode::qfloat<N>>( \
    const Matrix<newcode::qfloat<N>>&, const Params&, \
    std::vector<TileEarlyStopCounter>*, bool, const Matrix<float>*); \
template void process_window<newcode::qfloat<N>>( \
    Matrix<newcode::qfloat<N>>&, const Matrix<newcode::qfloat<N>>&, \
    std::size_t, std::size_t, const Params&, \
    std::size_t, std::size_t, std::size_t, \
    std::vector<TileEarlyStopCounter>*, bool, \
    const Matrix<float>*);

INSTANTIATE_QFLOAT(2)
INSTANTIATE_QFLOAT(3)
INSTANTIATE_QFLOAT(4)
INSTANTIATE_QFLOAT(5)
INSTANTIATE_QFLOAT(6)
INSTANTIATE_QFLOAT(7)
INSTANTIATE_QFLOAT(8)
INSTANTIATE_QFLOAT(9)
INSTANTIATE_QFLOAT(10)
INSTANTIATE_QFLOAT(11)
INSTANTIATE_QFLOAT(12)
INSTANTIATE_QFLOAT(13)
INSTANTIATE_QFLOAT(14)
INSTANTIATE_QFLOAT(15)

#undef INSTANTIATE_QFLOAT

template void process_window<float >(Matrix<float>&,  const Matrix<float>&,
                                     std::size_t, std::size_t, const Params&,
                                     std::size_t, std::size_t, std::size_t,
                                     std::vector<TileEarlyStopCounter>*, bool,
                                     const Matrix<float>*);
template void process_window<int8_t>(Matrix<int8_t>&, const Matrix<int8_t>&,
                                     std::size_t, std::size_t, const Params&,
                                     std::size_t, std::size_t, std::size_t,
                                     std::vector<TileEarlyStopCounter>*, bool,
                                     const Matrix<float>*);

} // namespace newcode
