#pragma once

#include "new_float_only/ofec_decoder.hpp"
#include "new_float_only/params.hpp"
#include "new_float_only/rx/ofec/chase/chase256.hpp" // 保留 Chase 头；本文档内有三参前向声明
#include "new_float_only/ofec_decoder_hard.hpp"
#include "new_float_only/rx/ofec/chase/decoder_core.hpp"

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

namespace new_float_only {
namespace detail {

using CoreFn = chase::DecoderCoreResult (*)(const matrix::Matrix<float>&,
                                            const matrix::Matrix<float>&,
                                            bool,
                                            const new_float_only::Params&);

} // namespace detail
} // namespace new_float_only

#include "ofec_tile_input.ipp"
#include "ofec_tile_decode.ipp"
#include "ofec_tile_writeback.ipp"

namespace new_float_only {
namespace detail {

/**
 * 单个 Tile 的完整处理流程：
 * 1. 重排成 256 维 Chase 输入；
 * 2. 调用 row core；
 * 3. 把 extrinsic 写回到 Tile 布局；
 * 4. 可选记录最后一个 soft Tile 的 history。
 */
TileProcessResult process_tile_impl(const matrix::Matrix<float>& tile_in,
                                    const matrix::Matrix<float>& ch_tile,
                                    const new_float_only::Params& p,
                                    size_t tile_top_row_global,
                                    bool use_hard_decode,
                                    bool normalize_extrinsic,
                                    const matrix::Matrix<float>* tx_llr_ref,
                                    CoreFn core_fn,
                                    matrix::Matrix<float>* last_tile_history_accum,
                                    bool capture_last_tile_history)
{
  constexpr int B         = static_cast<int>(new_float_only::Params::BITS_PER_SUBBLOCK_DIM);            // 16
  constexpr int N         = static_cast<int>(new_float_only::Params::NUM_SUBBLOCK_COLS * B);            // 128

  const size_t H = tile_in.rows();
  const size_t W = tile_in.cols();
  assert(W == static_cast<size_t>(N));
  assert(ch_tile.rows() == H && ch_tile.cols() == W);

  matrix::Matrix<float> tile_out = tile_in;

  const int SBR = p.CHASE_SBR;
  if (SBR != 1 && SBR != 2)
      throw std::invalid_argument("process_tile: CHASE_SBR must be 1 or 2.");

  const size_t rows_to_decode = static_cast<size_t>(SBR) * static_cast<size_t>(B);
  if (rows_to_decode == 0) {
      return TileProcessResult{tile_out};
  }

  TilePrepared prep = prepare_tile_inputs(tile_in, ch_tile, p,
                                          tile_top_row_global,
                                          SBR,
                                          rows_to_decode,
                                          tx_llr_ref);

  auto decoder_res = decode_tile(prep,
                                 use_hard_decode,
                                 normalize_extrinsic,
                                 p,
                                 core_fn);

  writeback_tile(prep,
                 decoder_res,
                 p,
                 tile_top_row_global,
                 capture_last_tile_history,
                 &tile_out,
                 last_tile_history_accum);

  return TileProcessResult{std::move(tile_out)};
}

} // namespace detail
} // namespace new_float_only
