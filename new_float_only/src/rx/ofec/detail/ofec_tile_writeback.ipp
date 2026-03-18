#pragma once

#include "ofec_tile_input.ipp"

#include <cassert>
#include <iostream>

namespace new_float_only {
namespace detail {

/**
 * 把 256 维行级 extrinsic 写回到 Tile 的 128 列布局中。
 * 写回规则与 prepare_tile_inputs 完全对偶：右半 128 位回到当前行，
 * 左半 128 位回到历史位置。若 capture_last_tile_history 为真，还会累计 history。
 */
void writeback_tile(const TilePrepared& prep,
                    const chase::DecoderCoreResult& decoder_res,
                    const new_float_only::Params& p,
                    size_t tile_top_row_global,
                    bool capture_last_tile_history,
                    matrix::Matrix<float>* tile_out,
                    matrix::Matrix<float>* last_tile_history_accum)
{
  constexpr int B         = static_cast<int>(new_float_only::Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(new_float_only::Params::NUM_SUBBLOCK_COLS * B);            // 128
  constexpr int K         = static_cast<int>(new_float_only::Params::BCH_K);                            // 239
  constexpr int TAKE_BITS = K - N;                                                      // 111
  constexpr int BCH_PAR   = static_cast<int>(new_float_only::Params::BCH_PARITY_BITS);                  // 16
  constexpr int OVR_IDX   = static_cast<int>(new_float_only::Params::BCH_OVERALL_IDX);                  // 255

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
      const int k = N + i;
      const size_t Ct = static_cast<size_t>((k - N) / B);
      const size_t ct = static_cast<size_t>((k % B) ^ r);
      const size_t col = Ct * static_cast<size_t>(B) + ct;
      const float extrinsic_llr = lout_row[static_cast<size_t>(k)];
      if (row_produced) {
        (*tile_out)[row_local][col] = extrinsic_llr;
      }
    }
    for (int j = 0; j < BCH_PAR; ++j) {
      const int k = K + j;
      const size_t Ct = static_cast<size_t>((k - N) / B);
      const size_t ct = static_cast<size_t>((k % B) ^ r);
      const size_t col = Ct * static_cast<size_t>(B) + ct;
      const float extrinsic_llr = lout_row[static_cast<size_t>(k)];
      if (row_produced) {
        (*tile_out)[row_local][col] = extrinsic_llr;
      }
    }
    {
      const int k = OVR_IDX;
      const size_t Ct = static_cast<size_t>((k - N) / B);
      const size_t ct = static_cast<size_t>((k % B) ^ r);
      const size_t col = Ct * static_cast<size_t>(B) + ct;
      const float extrinsic_llr = lout_row[static_cast<size_t>(k)];
      if (row_produced) {
        (*tile_out)[row_local][col] = extrinsic_llr;
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

      if (prep.trace.should_trace(prep.trace.trace_cfg.log_write_mapping, rr_global, cc_global)) {
        std::cout << " WRITE Mapping k=" << k
                  << " to global pos (" << rr_global << "," << cc_global << ")" << '\n';
        std::cout << "  Value=" << lout_row[static_cast<size_t>(k)] << '\n';
      }

      const long rr_local2 = rr_global - static_cast<long>(tile_top_row_global);
      const long cc_local2 = cc_global;

      const bool in_range =
          (rr_local2 >= 0 && rr_local2 < static_cast<long>(H) &&
           cc_local2 >= 0 && cc_local2 < static_cast<long>(W));

      assert(in_range && "process_tile: write-back out of tile range");

      const size_t rr_idx_local = static_cast<size_t>(rr_local2);
      const size_t cc_idx_local = static_cast<size_t>(cc_local2);
      const float prior_llr =
          (*tile_out)[rr_idx_local][cc_idx_local];
      const float extrinsic_llr = lout_row[static_cast<size_t>(k)];
      if (row_produced) {
        (*tile_out)[rr_idx_local][cc_idx_local] = extrinsic_llr;
      }

      if (capture_last_tile_history && last_tile_history_accum) {
        if (rr_global >= 0 && cc_global >= 0) {
          const size_t rr_idx_global = static_cast<size_t>(rr_global);
          const size_t cc_idx_global = static_cast<size_t>(cc_global);
          if (rr_idx_global < last_tile_history_accum->rows() &&
              cc_idx_global < last_tile_history_accum->cols()) {
            // 这里保留旧 float 路径的 history 口径：新 extrinsic 与旧 history 直接相加。
            (*last_tile_history_accum)[rr_idx_global][cc_idx_global] =
                extrinsic_llr + prior_llr;
          }
        }
      }
    }
  }
}

} // namespace detail
} // namespace new_float_only
