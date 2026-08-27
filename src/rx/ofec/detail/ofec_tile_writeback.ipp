#pragma once

#include "ofec_tile_input.ipp"

namespace newcode {
namespace detail {

template <typename LLR>
void writeback_tile(const TilePrepared<LLR>& prep,
                    const chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>& decoder_res,
                    const newcode::Params& p,
                    size_t tile_top_row_global,
                    bool capture_last_tile_history,
                    matrix::Matrix<LLR>* tile_out,
                    matrix::Matrix<float>* last_tile_history_accum,
                    bool preserve_history_when_not_produced = true)
{
  // 输入:
  // - prep: tile 输入准备阶段留下的映射关系和 trace 信息。
  // - decoder_res: decoder core 输出，包括 lout 和 produced_rows。
  // - p: 当前 tile 参数。
  // - tile_top_row_global: tile 顶部的全局行号。
  // - capture_last_tile_history: 是否记录当前 tile 的历史值。
  // - tile_out: tile 输出矩阵，会被原地写回。
  // - last_tile_history_accum: 可选的全局历史矩阵。
  // 输出:
  // - 无返回值；通过 tile_out / last_tile_history_accum 反映写回结果。
  // 用途:
  // - 把 decoder core 的线性输出重新映射回 tile 的二维布局，并可选累积历史值。
  using Adapter = LinMatrixAdapter<LLR>;
  using CoreLLR = typename Adapter::core_type;

  constexpr int B         = static_cast<int>(newcode::Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(newcode::Params::NUM_SUBBLOCK_COLS * B);            // 128
  constexpr int K         = static_cast<int>(newcode::Params::BCH_K);                            // 239
  constexpr int TAKE_BITS = K - N;                                                      // 111
  constexpr int BCH_PAR   = static_cast<int>(newcode::Params::BCH_PARITY_BITS);                  // 16
  constexpr int OVR_IDX   = static_cast<int>(newcode::Params::BCH_OVERALL_IDX);                  // 255

  auto core_to_float = [](const CoreLLR& value) -> float {
    // 将 core 域类型统一转成 float，便于做历史累积。
    if constexpr (std::is_same_v<CoreLLR, float> || std::is_same_v<CoreLLR, double>) {
      return static_cast<float>(value);
    } else {
      return llr_to_float(value);
    }
  };
  constexpr bool history_needs_dequant =
      std::is_same_v<CoreLLR, float> &&
      !std::is_floating_point_v<LLR> &&
      !std::is_integral_v<LLR>;
  auto history_value = [&](const CoreLLR& combined,
                           float extrinsic,
                           const LLR& prior) -> float {
    // 对 qfloat 路径，combined 可能仍在 code 域，因此历史值需要反量化后再保存。
    if constexpr (history_needs_dequant) {
      // For quantized LLR (e.g., qfloat::qfloat<N>), combined is in code domain.
      // Sum dequantized values so history stores real amplitudes.
      return extrinsic + llr_to_float(prior);
    } else {
      return core_to_float(combined);
    }
  };

  const size_t H = tile_out->rows();
  const size_t W = tile_out->cols();

  for (size_t row_idx = 0; row_idx < prep.row_local_lookup.size(); ++row_idx)
  {
    if (row_idx >= decoder_res.produced_rows.size()) continue;
    const bool row_produced = decoder_res.produced_rows[row_idx];

    const size_t row_local  = prep.row_local_lookup[row_idx];
    const size_t row_global = prep.row_global_lookup[row_idx];
    const auto&  lout_row   = decoder_res.lout[row_idx];

    const int r = static_cast<int>(row_global % static_cast<size_t>(B));

    for (int i = 0; i < TAKE_BITS; ++i) {
      // 先写回当前行上的系统/新信息位。
      const int k = N + i;
      const size_t Ct = static_cast<size_t>((k - N) / B);
      const size_t ct = static_cast<size_t>((k % B) ^ r);
      const size_t col = Ct * static_cast<size_t>(B) + ct;
      const LLR prior_llr = (*tile_out)[row_local][col];
      const LLR extrinsic_llr =
          qfloat::llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      if (row_produced) {
        (*tile_out)[row_local][col] = extrinsic_llr;
      }


    }
    for (int j = 0; j < BCH_PAR; ++j) {
      // 再写回 BCH parity。
      const int k = K + j;
      const size_t Ct = static_cast<size_t>((k - N) / B);
      const size_t ct = static_cast<size_t>((k % B) ^ r);
      const size_t col = Ct * static_cast<size_t>(B) + ct;
      const LLR prior_llr = (*tile_out)[row_local][col];
      const LLR extrinsic_llr =
           qfloat::llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      if (row_produced) {
        (*tile_out)[row_local][col] = extrinsic_llr;
      }

    }
    {
      // 最后写回 overall parity。
      const int k = OVR_IDX;
      const size_t Ct = static_cast<size_t>((k - N) / B);
      const size_t ct = static_cast<size_t>((k % B) ^ r);
      const size_t col = Ct * static_cast<size_t>(B) + ct;
      const LLR prior_llr = (*tile_out)[row_local][col];
      const LLR extrinsic_llr =
          qfloat::llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      if (row_produced) {
        (*tile_out)[row_local][col] = extrinsic_llr;
      }


    }

    const long R = static_cast<long>(row_global / static_cast<size_t>(B));
    for (int k = 0; k < N; ++k)
    {
      // 前 N 位需要重新投影回“历史信息”所在的全局二维坐标。
      const long br = (R ^ 1L)
                    - static_cast<long>(2 * p.NUM_GUARD_SUBROWS)
                    - static_cast<long>(2 * (N / B))
                    + static_cast<long>(2 * (k / B));
      const long bc = static_cast<long>(k / B);
      const long bit_row_in_block = static_cast<long>((k % B) ^ r);
      const long bit_col_in_block = static_cast<long>(r);

      const long rr_global = br * B + bit_row_in_block;
      const long cc_global = bc * B + bit_col_in_block;

      if (prep.trace.should_trace(prep.trace.trace_cfg.log_write_mapping, rr_global, cc_global)) {
        std::cout << " WRITE Mapping k=" << k
                  << " to global pos (" << rr_global << "," << cc_global << ")" << '\n';
        std::cout << "  Value=" << qfloat::llr_to_float(lout_row[static_cast<size_t>(k)]) << '\n';
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
          (*tile_out)[rr_idx_local][cc_idx_local];
      LLR extrinsic_llr =
          qfloat::llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      if (row_produced) {
        (*tile_out)[rr_idx_local][cc_idx_local] = extrinsic_llr;
      }

      if (capture_last_tile_history && last_tile_history_accum) {
        if (rr_global >= 0 && cc_global >= 0) {
          const size_t rr_idx_global = static_cast<size_t>(rr_global);
          const size_t cc_idx_global = static_cast<size_t>(cc_global);
          if (rr_idx_global < last_tile_history_accum->rows() &&
              cc_idx_global < last_tile_history_accum->cols()) {
            if (row_produced) {
              const float extrinsic_llr_last_tile =
                  lout_row[static_cast<size_t>(k)] / p.ALPHA;
              // 已调度行保存“新 extrinsic + 旧 prior”。
              const CoreLLR combined =
                  Adapter::combine(extrinsic_llr, prior_llr);
              (*last_tile_history_accum)[rr_idx_global][cc_idx_global] =
                  history_value(combined, extrinsic_llr_last_tile, prior_llr);
            } else if (preserve_history_when_not_produced) {
              // 未调度行没有新 extrinsic，但最终 history 必须透传当前 prior。
              // Buffered FIFO 的 pending/AlreadyDecoded 行在后续服务时刻没有
              // 产生新结果；此时由调用方关闭该透传，避免覆盖同一物理行已经
              // 写入的最终 Level 6 history。
              (*last_tile_history_accum)[rr_idx_global][cc_idx_global] =
                  qfloat::llr_to_float(prior_llr);
            }
          }
        }
      }
    }
  }
}

} // namespace detail
} // namespace newcode
