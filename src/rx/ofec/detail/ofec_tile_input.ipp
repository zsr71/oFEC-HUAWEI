#pragma once

namespace newcode {
namespace detail {

struct TileTraceContext {
  newcode::Params::DebugTraceConfig trace_cfg;
  bool trace_has_coords{false};
  bool trace_has_targets{false};
  bool trace_enabled{false};
  long trace_row{-1};
  long trace_col{-1};

  bool should_trace(bool flag, long rr, long cc) const {
    // 输入:
    // - flag: 当前类别日志是否开启。
    // - rr/cc: 待判断的全局坐标。
    // 输出:
    // - 是否需要对这个坐标打印/记录 trace。
    // 用途:
    // - 集中封装 trace 条件判断，避免散落在映射代码里。
    return trace_enabled && flag && rr == trace_row && cc == trace_col;
  }
};

template <typename LLR>
struct TilePrepared {
  using Adapter = LinMatrixAdapter<LLR>;
  using CoreLLR = typename Adapter::core_type;

  matrix::Matrix<CoreLLR> lin_matrix;
  matrix::Matrix<CoreLLR> lch_matrix;
  std::vector<size_t> row_local_lookup;
  std::vector<size_t> row_global_lookup;
  newcode::Params params_for_core;
  std::shared_ptr<std::vector<std::vector<int8_t>>> expected_bits;
  TileTraceContext trace;
};

template <typename LLR>
TilePrepared<LLR> prepare_tile_inputs(const matrix::Matrix<LLR>& tile_in,
                                      const matrix::Matrix<LLR>& ch_tile,
                                      const newcode::Params& p,
                                      size_t tile_top_row_global,
                                      int SBR,
                                      size_t rows_to_decode,
                                      const matrix::Matrix<float>* tx_llr_ref)
{
  // 输入:
  // - tile_in: 当前 tile 的先验/历史矩阵。
  // - ch_tile: 当前 tile 的信道矩阵。
  // - p: 当前 tile 参数。
  // - tile_top_row_global: tile 顶部全局行号。
  // - SBR: 本 tile 需要解码的 subblock row 数。
  // - rows_to_decode: 实际 decoder row 数，通常为 SBR * B。
  // - tx_llr_ref: 可选参考矩阵，仅用于调试时生成 expected bit。
  // 输出:
  // - TilePrepared，包含 lin_matrix/lch_matrix、行映射表、调试信息和 core 参数。
  // 用途:
  // - 将 tile 的二维几何布局重排成 decoder core 需要的 2N 线性输入格式。
  using Adapter = typename TilePrepared<LLR>::Adapter;
  constexpr int B         = static_cast<int>(newcode::Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(newcode::Params::NUM_SUBBLOCK_COLS * B);            // 128
  constexpr int K         = static_cast<int>(newcode::Params::BCH_K);                            // 239
  constexpr int TAKE_BITS = K - N;                                                      // 111
  constexpr int BCH_PAR   = static_cast<int>(newcode::Params::BCH_PARITY_BITS);                  // 16
  constexpr int OVR_IDX   = static_cast<int>(newcode::Params::BCH_OVERALL_IDX);                  // 255

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

  newcode::Params params_for_core = p;
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

  matrix::Matrix<typename LinMatrixAdapter<LLR>::core_type> lin_matrix(rows_to_decode,
                                                               static_cast<size_t>(2 * N));
  matrix::Matrix<typename LinMatrixAdapter<LLR>::core_type> lch_matrix(rows_to_decode,
                                                               static_cast<size_t>(2 * N));
  std::vector<size_t> row_local_lookup(rows_to_decode, 0);
  std::vector<size_t> row_global_lookup(rows_to_decode, 0);

  auto should_trace = [&](bool flag, long rr, long cc) -> bool {
    return trace_ctx.should_trace(flag, rr, cc);
  };

  auto register_chase_entry = [&](const std::string& label,
                                  long bit_index,
                                  long rr, long cc,
                                  size_t row_idx, int k) {
    // 为被追踪的全局坐标建立“decoder 行/列 <-> 全局位置”的映射记录。
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

  auto update_chase_entry_channel_llr = [&](long rr, long cc,
                                            size_t row_idx, int k,
                                            const LLR& lch) {
    auto& entries = params_for_core.debug_trace.active_chase_entries;
    for (auto& entry : entries) {
      if (entry.row_index != static_cast<int>(row_idx)) continue;
      if (entry.k != k) continue;
      if (entry.global_row != rr || entry.global_col != cc) continue;
      entry.has_channel_llr_float = true;
      entry.channel_llr_float = qfloat::llr_to_float(lch);
      if constexpr (!std::is_floating_point_v<LLR> && !std::is_integral_v<LLR>) {
        entry.has_channel_llr_code = true;
        entry.channel_llr_code = lch.code();
      } else {
        entry.has_channel_llr_code = false;
      }
    }
  };

  for (int s = 0; s < SBR; ++s){
    // s 表示当前处理的第几个 subblock row 组。
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
        // 前 N 位对应“历史信息”部分，要从上一轮相关位置重新取值并做交织映射。
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
          const LLR Lch = ch_tile[static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)];
          const LLR La  = tile_in [static_cast<size_t>(rr_local2)][static_cast<size_t>(cc_local2)];
          // lin_matrix 是送给 core 的总输入；lch_matrix 只保留信道项。
          lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
          lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
          update_chase_entry_channel_llr(rr_global, cc_global, row_idx, k, Lch);
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
        // 接下来的 TAKE_BITS 部分对应当前行上的系统/新信息位。
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
        const LLR Lch = ch_tile[row_local][src_col];
        const LLR La  = tile_in [row_local][src_col];
        lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
        lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
        update_chase_entry_channel_llr(static_cast<long>(row_local + tile_top_row_global),
                                       static_cast<long>(src_col), row_idx, k, Lch);
        if (expected_bits &&
            (tile_top_row_global + row_local) < tx_llr_ref->rows() &&
            src_col < tx_llr_ref->cols()) {
          const float v = (*tx_llr_ref)[tile_top_row_global + row_local][src_col];
          (*expected_bits)[row_idx][static_cast<size_t>(k)] = (v < 0.0f) ? 1 : 0;
        }
      }

      for (int j = 0; j < BCH_PAR; ++j) {
        // 然后是 BCH 奇偶校验位。
        const int k = K + j;
        const size_t Ct = static_cast<size_t>((k - N) / B);
        const size_t ct = static_cast<size_t>((k % B) ^ r);
        const size_t src_col = Ct * static_cast<size_t>(B) + ct;

        const LLR Lch = ch_tile[row_local][src_col];
        const LLR La  = tile_in [row_local][src_col];
        lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
        lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
        update_chase_entry_channel_llr(static_cast<long>(row_local + tile_top_row_global),
                                       static_cast<long>(src_col), row_idx, k, Lch);
        if (expected_bits &&
            (tile_top_row_global + row_local) < tx_llr_ref->rows() &&
            src_col < tx_llr_ref->cols()) {
          const float v = (*tx_llr_ref)[tile_top_row_global + row_local][src_col];
          (*expected_bits)[row_idx][static_cast<size_t>(k)] = (v < 0.0f) ? 1 : 0;
        }
      }

      {
        // 最后一位是 overall parity。
        const int k = OVR_IDX;
        const size_t Ct = static_cast<size_t>((k - N) / B);
        const size_t ct = static_cast<size_t>((k % B) ^ r);
        const size_t src_col = Ct * static_cast<size_t>(B) + ct;

        const LLR Lch = ch_tile[row_local][src_col];
        const LLR La  = tile_in [row_local][src_col];
        lin_matrix[row_idx][static_cast<size_t>(k)] = Adapter::combine(Lch, La);
        lch_matrix[row_idx][static_cast<size_t>(k)] = Adapter::channel(Lch);
        update_chase_entry_channel_llr(static_cast<long>(row_local + tile_top_row_global),
                                       static_cast<long>(src_col), row_idx, k, Lch);
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
    // register_chase_entry 可能发生在 expected_bits 该位置真正写入之前，
    // 这里在整块 lin_matrix 构造完成后，统一把活动追踪项的 expected_bit 回填成最终值。
    for (auto& entry : params_for_core.debug_trace.active_chase_entries) {
      if (entry.row_index < 0 || entry.k < 0) continue;
      const size_t row_idx = static_cast<size_t>(entry.row_index);
      const size_t k_idx = static_cast<size_t>(entry.k);
      if (row_idx >= expected_bits->size()) continue;
      if (k_idx >= (*expected_bits)[row_idx].size()) continue;
      entry.expected_bit = (*expected_bits)[row_idx][k_idx];
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

  TilePrepared<LLR> prep;
  // 将所有为 core 准备好的中间结果打包返回。
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
} // namespace newcode
