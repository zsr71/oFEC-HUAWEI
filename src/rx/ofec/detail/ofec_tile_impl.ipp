#pragma once

#include "newcode/ofec_decoder.hpp"
#include "newcode/params.hpp"
#include "newcode/rx/ofec/chase/chase256.hpp" // 保留 Chase 头；本文档内有三参前向声明
#include "newcode/ofec_decoder_hard.hpp"
#include "newcode/rx/ofec/chase/decoder_core.hpp"
#include "newcode/llr_utils.hpp"
#include "newcode/decoder_api.hpp"
#include "newcode/ofec/common/lin_matrix_adapters.hpp"
#include "newcode/ofec/earltstop/tile_early_stop_stats.hpp"

#include <filesystem>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cctype>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace newcode {
namespace detail {

template <typename LLR>
using CoreFn = chase::DecoderCoreResult<LLR> (*)(const matrix::Matrix<LLR>&,
                                                 const matrix::Matrix<LLR>&,
                                                 bool,
                                                 const newcode::Params&,
                                                 const std::vector<bool>* early_stop_row_flags);

} // namespace detail
} // namespace newcode

#include "ofec_tile_input.ipp"
#include "ofec_tile_decode.ipp"
#include "ofec_tile_writeback.ipp"

namespace newcode {
namespace detail {

template <typename LLR>
TileProcessResult<LLR> process_tile_impl(const matrix::Matrix<LLR>& tile_in,
                                         const matrix::Matrix<LLR>& ch_tile,
                                         const newcode::Params& p,
                                         size_t tile_top_row_global,
                                         bool use_hard_decode,
                                         bool normalize_extrinsic,
                                         const matrix::Matrix<float>* tx_llr_ref,
                                         CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
                                         matrix::Matrix<float>* last_tile_history_accum,
                                         bool capture_last_tile_history)
{
  constexpr int B         = static_cast<int>(newcode::Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(newcode::Params::NUM_SUBBLOCK_COLS * B);            // 128

  const size_t H = tile_in.rows();
  const size_t W = tile_in.cols();
  assert(W == static_cast<size_t>(N));
  assert(ch_tile.rows() == H && ch_tile.cols() == W);

  matrix::Matrix<LLR> tile_out = tile_in;

  const int SBR = p.CHASE_SBR;
  if (SBR != 1 && SBR != 2)
      throw std::invalid_argument("process_tile: CHASE_SBR must be 1 or 2.");

  const size_t rows_to_decode = static_cast<size_t>(SBR) * static_cast<size_t>(B);
  if (rows_to_decode == 0) {
      return TileProcessResult<LLR>{tile_out, false, 0, 0};
  }

  TilePrepared<LLR> prep = prepare_tile_inputs(tile_in, ch_tile, p,
                                               tile_top_row_global,
                                               SBR,
                                               rows_to_decode,
                                               tx_llr_ref);

  TileEarlyStopResult early_stop_stats = tile_early_stop_stats(prep.lin_matrix);
  bool early_stop_triggered = early_stop_stats.all_rows_passed;

  auto decoder_res = decode_tile<LLR>(prep,
                                      use_hard_decode,
                                      normalize_extrinsic,
                                      p,
                                      &early_stop_stats.row_passed_flags,
                                      core_fn);

  writeback_tile(prep,
                 decoder_res,
                 p,
                 tile_top_row_global,
                 capture_last_tile_history,
                 &tile_out,
                 last_tile_history_accum);

  return TileProcessResult<LLR>{
      std::move(tile_out),
      early_stop_triggered,
      early_stop_stats.rows_passed,
      early_stop_stats.rows_total};
}

} // namespace detail
} // namespace newcode
