#pragma once

#include "ofec_tile_input.ipp"

namespace newcode {
namespace detail {

template <typename LLR>
void writeback_tile(const TilePrepared<LLR>& prep,
                    const DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type>& decoder_res,
                    const Params& p,
                    size_t tile_top_row_global,
                    bool capture_last_tile_history,
                    Matrix<LLR>* tile_out,
                    Matrix<float>* last_tile_history_accum)
{
  using Adapter = LinMatrixAdapter<LLR>;
  using CoreLLR = typename Adapter::core_type;

  constexpr int B         = static_cast<int>(Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(Params::NUM_SUBBLOCK_COLS * B);            // 128
  constexpr int K         = static_cast<int>(Params::BCH_K);                            // 239
  constexpr int TAKE_BITS = K - N;                                                      // 111
  constexpr int BCH_PAR   = static_cast<int>(Params::BCH_PARITY_BITS);                  // 16
  constexpr int OVR_IDX   = static_cast<int>(Params::BCH_OVERALL_IDX);                  // 255

  auto core_to_float = [](const CoreLLR& value) -> float {
    if constexpr (std::is_same_v<CoreLLR, float> || std::is_same_v<CoreLLR, double>) {
      return static_cast<float>(value);
    } else {
      return llr_to_float(value);
    }
  };

  const size_t H = tile_out->rows();
  const size_t W = tile_out->cols();

  for (size_t row_idx = 0; row_idx < prep.row_local_lookup.size(); ++row_idx)
  {
    if (row_idx >= decoder_res.produced_rows.size()) continue;
    if (!decoder_res.produced_rows[row_idx]) continue;

    const size_t row_local  = prep.row_local_lookup[row_idx];
    const size_t row_global = prep.row_global_lookup[row_idx];
    const auto&  lout_row   = decoder_res.lout[row_idx];

    const int r = static_cast<int>(row_global % static_cast<size_t>(B));

    for (int i = 0; i < TAKE_BITS; ++i) {
      const int k = N + i;
      const size_t Ct = static_cast<size_t>((k - N) / B);
      const size_t ct = static_cast<size_t>((k % B) ^ r);
      const size_t col = Ct * static_cast<size_t>(B) + ct;
      const LLR prior_llr = (*tile_out)[row_local][col];
      const LLR extrinsic_llr =
          llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      (*tile_out)[row_local][col] = extrinsic_llr;

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
      const LLR prior_llr = (*tile_out)[row_local][col];
      const LLR extrinsic_llr =
          llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      (*tile_out)[row_local][col] = extrinsic_llr;

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
      const LLR prior_llr = (*tile_out)[row_local][col];
      const LLR extrinsic_llr =
          llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      (*tile_out)[row_local][col] = extrinsic_llr;

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

      if (prep.trace.should_trace(prep.trace.trace_cfg.log_write_mapping, rr_global, cc_global)) {
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
          (*tile_out)[rr_idx_local][cc_idx_local];
      const LLR extrinsic_llr =
          llr_from_float<LLR>(lout_row[static_cast<size_t>(k)]);
      (*tile_out)[rr_idx_local][cc_idx_local] = extrinsic_llr;

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

} // namespace detail
} // namespace newcode
