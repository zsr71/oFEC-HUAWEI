#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace new_float_only {
namespace detail {

struct TileTraceContext {
  new_float_only::Params::DebugTraceConfig trace_cfg;
  bool trace_has_coords{false};
  bool trace_has_targets{false};
  bool trace_enabled{false};
  long trace_row{-1};
  long trace_col{-1};

  bool should_trace(bool flag, long rr, long cc) const {
    return trace_enabled && flag && rr == trace_row && cc == trace_col;
  }
};

struct TilePrepared {
  matrix::Matrix<float> lin_matrix;
  matrix::Matrix<float> lch_matrix;
  std::vector<size_t> row_local_lookup;
  std::vector<size_t> row_global_lookup;
  new_float_only::Params params_for_core;
  std::shared_ptr<std::vector<std::vector<int8_t>>> expected_bits;
  TileTraceContext trace;
};

/**
 * 把当前 Tile 的 128 列工作矩阵重排为 Chase 所需的 256 维行码字输入。
 * 关键映射：
 * 1. k<128 时从历史区按 oFEC 历史公式读取；
 * 2. k>=128 时从当前 Tile 底部待译码行直接抽取系统位/校验位/overall parity。
 */
TilePrepared prepare_tile_inputs(const matrix::Matrix<float>& tile_in,
                                 const matrix::Matrix<float>& ch_tile,
                                 const new_float_only::Params& p,
                                 size_t tile_top_row_global,
                                 int SBR,
                                 size_t rows_to_decode,
                                 const matrix::Matrix<float>* tx_llr_ref)
{
  constexpr int B         = static_cast<int>(new_float_only::Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(new_float_only::Params::NUM_SUBBLOCK_COLS * B);            // 128
  constexpr int K         = static_cast<int>(new_float_only::Params::BCH_K);                            // 239
  constexpr int TAKE_BITS = K - N;                                                      // 111
  constexpr int BCH_PAR   = static_cast<int>(new_float_only::Params::BCH_PARITY_BITS);                  // 16
  constexpr int OVR_IDX   = static_cast<int>(new_float_only::Params::BCH_OVERALL_IDX);                  // 255

  const size_t H = tile_in.rows();
  const size_t W = tile_in.cols();
  (void)W;

  TileTraceContext trace_ctx;
  trace_ctx.trace_cfg = p.debug_trace;
  trace_ctx.trace_has_coords = (trace_ctx.trace_cfg.row >= 0 && trace_ctx.trace_cfg.col >= 0);
  trace_ctx.trace_has_targets = !trace_ctx.trace_cfg.targets.empty();
  trace_ctx.trace_enabled =
      trace_ctx.trace_cfg.enable &&
      (trace_ctx.trace_has_coords || trace_ctx.trace_has_targets);
  trace_ctx.trace_row = trace_ctx.trace_has_coords ? trace_ctx.trace_cfg.row : -1;
  trace_ctx.trace_col = trace_ctx.trace_has_coords ? trace_ctx.trace_cfg.col : -1;

  new_float_only::Params params_for_core = p;
  params_for_core.debug_trace.active_chase_entries.clear();
  if (params_for_core.debug_trace.chase_expected_bits) {
    params_for_core.debug_trace.chase_expected_bits.reset();
    params_for_core.debug_trace.chase_expected_bits_row = nullptr;
  }

  std::shared_ptr<std::vector<std::vector<int8_t>>> expected_bits;
  if (tx_llr_ref) {
    expected_bits = std::make_shared<std::vector<std::vector<int8_t>>>(
        rows_to_decode, std::vector<int8_t>(static_cast<size_t>(2 * N), -1));
  }

  matrix::Matrix<float> lin_matrix(rows_to_decode, static_cast<size_t>(2 * N));
  matrix::Matrix<float> lch_matrix(rows_to_decode, static_cast<size_t>(2 * N));
  std::vector<size_t> row_local_lookup(rows_to_decode, 0);
  std::vector<size_t> row_global_lookup(rows_to_decode, 0);

  auto should_trace = [&](bool flag, long rr, long cc) -> bool {
    return trace_ctx.should_trace(flag, rr, cc);
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
    new_float_only::Params::DebugTraceConfig::ChaseTraceEntry entry;
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
    if (!trace_ctx.trace_enabled) return;
    if (trace_ctx.trace_has_coords && rr == trace_ctx.trace_row && cc == trace_ctx.trace_col) {
      const std::string label =
          "row" + std::to_string(rr) + "_col" + std::to_string(cc);
      register_chase_entry(label, -1, rr, cc, row_idx, k);
    }
    for (size_t idx = 0; idx < trace_ctx.trace_cfg.targets.size(); ++idx) {
      const auto& target = trace_ctx.trace_cfg.targets[idx];
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

        if (should_trace(trace_ctx.trace_cfg.log_read_mapping, rr_global, cc_global)) {
          std::cout << " AS History READ Mapping k=" << k
                    << " to global pos (" << rr_global << "," << cc_global << ")" << '\n';
        }
        try_mark_chase_coord(rr_global, cc_global, row_idx, k);

        const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
        const long cc_local2 = cc_global;

        if (rr_local2 >= 0 && rr_local2 < static_cast<long>(H) &&
            cc_local2 >= 0 && cc_local2 < static_cast<long>(W))
        {
          const float Lch = ch_tile[static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)];
          const float La  = tile_in [static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)];
          lin_matrix[row_idx][static_cast<size_t>(k)] = Lch + La;
          lch_matrix[row_idx][static_cast<size_t>(k)] = Lch;
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
        if (should_trace(trace_ctx.trace_cfg.log_read_mapping, row_local+tile_top_row_global, src_col)) {
          std::cout << " AS New Information READ Mapping k=" << k
                    << " to global pos (" << (row_local + tile_top_row_global) << "," << src_col << ")" << '\n';
        }
        try_mark_chase_coord(static_cast<long>(row_local + tile_top_row_global),
                             static_cast<long>(src_col), row_idx, k);
        const float Lch = ch_tile[row_local][src_col];
        const float La  = tile_in [row_local][src_col];
        lin_matrix[row_idx][static_cast<size_t>(k)] = Lch + La;
        lch_matrix[row_idx][static_cast<size_t>(k)] = Lch;
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

        const float Lch = ch_tile[row_local][src_col];
        const float La  = tile_in [row_local][src_col];
        lin_matrix[row_idx][static_cast<size_t>(k)] = Lch + La;
        lch_matrix[row_idx][static_cast<size_t>(k)] = Lch;
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

        const float Lch = ch_tile[row_local][src_col];
        const float La  = tile_in [row_local][src_col];
        lin_matrix[row_idx][static_cast<size_t>(k)] = Lch + La;
        lch_matrix[row_idx][static_cast<size_t>(k)] = Lch;
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

  TilePrepared prep;
  prep.lin_matrix = std::move(lin_matrix);
  prep.lch_matrix = std::move(lch_matrix);
  prep.row_local_lookup = std::move(row_local_lookup);
  prep.row_global_lookup = std::move(row_global_lookup);
  prep.params_for_core = std::move(params_for_core);
  prep.expected_bits = std::move(expected_bits);
  prep.trace = std::move(trace_ctx);
  return prep;
}

} // namespace detail
} // namespace new_float_only
