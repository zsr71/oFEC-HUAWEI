#pragma once

#include "ofec_window_impl.ipp"

#include "newcode/decoder_api.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/ofec/mux/mux_config_validate.hpp"
#include "newcode/ofec/mux/mux_group_config_validate.hpp"
#include "newcode/quantized_llr_dump.hpp"

#include <string>
#include <cmath>
#include <utility>
#include <vector>

namespace newcode {
namespace detail {

template <typename LLR>
matrix::Matrix<LLR> ofec_decode_llr_impl(const matrix::Matrix<LLR>& llr_mat, const newcode::Params& p,
                                 std::vector<TileEarlyStopCounter>* tile_stats,
                                 bool normalize_extrinsic,
                                 const matrix::Matrix<float>* tx_llr_ref,
                                 CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn)
{
  // 输入:
  // - llr_mat: 整帧输入 LLR 矩阵，按行组织。
  // - p: OFEC 解码参数，包含窗口、tile、alpha/beta、mux 等配置。
  // - tile_stats: 可选输出，用于累计每个 tile 的 early-stop 统计。
  // - normalize_extrinsic: 是否对 core 输出的外信息做归一化。
  // - tx_llr_ref: 可选参考矩阵，仅用于调试/追踪比特。
  // - core_fn: 真正的 Chase/Decoder core 回调。
  // 输出:
  // - 返回一张与输入同尺寸的 LLR 矩阵，表示窗口滑动和 tile 累积后的最终结果。
  // 用途:
  // - 这是 OFEC 软译码 detail 层的总入口，负责参数校验、窗口调度、
  //   工作矩阵初始化、窗口循环处理，以及最终把历史信息重新与信道项合成。
  const size_t N = newcode::Params::NUM_SUBBLOCK_COLS * newcode::Params::BITS_PER_SUBBLOCK_DIM;

  const size_t RROWS = llr_mat.rows();
  const size_t CCOLS = llr_mat.cols();
  if (CCOLS != N)
      throw std::invalid_argument("ofec_decode_llr: llr_mat cols != N.");

  assert(p.valid());
  if (!std::isfinite(p.TWOMAIN_HISO_M2) || p.TWOMAIN_HISO_M2 < 0.0f ||
      !std::isfinite(p.TWOMAIN_HISO_RHO_CORR) ||
      p.TWOMAIN_HISO_RHO_CORR < 0.0f ||
      p.TWOMAIN_HISO_RHO_CORR > 1.0f ||
      !std::isfinite(p.TWOMAIN_HISO_RHO_KEEP) ||
      p.TWOMAIN_HISO_RHO_KEEP < 0.0f ||
      p.TWOMAIN_HISO_RHO_KEEP > 1.0f) {
    throw std::invalid_argument(
        "ofec_decode_llr: invalid TwoMain HISO output parameters");
  }
  validate_level56_shared_config(p);
  const std::size_t mux_tile_count =
      p.LEVEL56_SHARED_ENABLE ? kLevel5TileIndex : p.TILES_PER_WIN;
  const auto mux_ok = p.LEVEL56_SHARED_ENABLE
      ? newcode::mux::validate_siso_active_prefix(p.SISO_ACTIVE_LIST,
                                                  mux_tile_count)
      : newcode::mux::validate_siso_active_list(p.SISO_ACTIVE_LIST,
                                                mux_tile_count);
  if (!mux_ok.ok) {
    throw std::invalid_argument("ofec_decode_llr: " + mux_ok.error);
  }
  const auto hiho_ok = p.LEVEL56_SHARED_ENABLE
      ? newcode::mux::validate_hiho_active_prefix(p.HIHO_ACTIVE_LIST,
                                                  mux_tile_count)
      : newcode::mux::validate_hiho_active_list(p.HIHO_ACTIVE_LIST,
                                                mux_tile_count);
  if (!hiho_ok.ok) {
    throw std::invalid_argument("ofec_decode_llr: " + hiho_ok.error);
  }
  if (!p.HYBRID_ENABLE_LIST.empty() &&
      p.HYBRID_ENABLE_LIST.size() != p.TILES_PER_WIN) {
    throw std::invalid_argument(
        "ofec_decode_llr: HYBRID_ENABLE_LIST length must equal TILES_PER_WIN");
  }
  if (!p.HYBRID_HARD_LLR_MAG_LIST.empty() &&
      p.HYBRID_HARD_LLR_MAG_LIST.size() != p.TILES_PER_WIN) {
    throw std::invalid_argument(
        "ofec_decode_llr: HYBRID_HARD_LLR_MAG_LIST length must equal TILES_PER_WIN");
  }
  const std::size_t rows_to_decode =
      static_cast<std::size_t>(p.CHASE_SBR) *
      newcode::Params::BITS_PER_SUBBLOCK_DIM;
  const auto group_ok =
      newcode::mux::validate_mux_group_g(p.MUX_GROUP_G, rows_to_decode);
  if (!group_ok.ok) {
    throw std::invalid_argument("ofec_decode_llr: " + group_ok.error);
  }
  if (p.MUX_ENABLE_RECONFIG) {
    for (std::size_t tile_idx = 0; tile_idx < mux_tile_count; ++tile_idx) {
      const auto reconfig_ok = newcode::mux::validate_mux_reconfig_runtime(
          p.MUX_GROUP_G, p.SISO_ACTIVE_LIST[tile_idx], rows_to_decode);
      if (!reconfig_ok.ok) {
        throw std::invalid_argument(
            "ofec_decode_llr: tile " + std::to_string(tile_idx) + ": " +
            reconfig_ok.error);
      }
    }
  }

  const size_t TILE_HEIGHT_ROWS = p.tile_height_rows();
  const size_t TILE_STRIDE_ROWS = p.tile_stride_rows();
  const size_t WIN_HEIGHT_ROWS  = p.win_height_rows();
  const size_t POP_PUSH_ROWS    = p.pop_push_rows();
  const size_t TILES_PER_WIN    = p.TILES_PER_WIN;


  matrix::Matrix<LLR> channel_llr = llr_mat;
  matrix::Matrix<LLR> work_llr(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
          // work_llr 保存“当前已累积的外信息”，初始为 0。
          work_llr[r][c] = qfloat::llr_from_float<LLR>(0.0f);
  matrix::Matrix<float> last_tile_history_llr(RROWS, N);

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
  std::size_t level56_shared_invocation = 0;
  Level56TemporalState<LLR> level56_temporal_state;
  Level56BufferedFifoState<LLR> level56_buffered_state;

  while (win_start <= last_ws) {
    const size_t win_end = win_start + WIN_HEIGHT_ROWS - 1;
    // 逐窗口处理；每个窗口内部会再拆成多个 tile。
    process_window_impl<LLR>(work_llr, channel_llr,
                             win_start, win_end, p,
                             TILE_HEIGHT_ROWS, TILE_STRIDE_ROWS, TILES_PER_WIN,
                             stats_ptr,
                             normalize_extrinsic,
                             tx_llr_ref,
                             core_fn,
                             &last_tile_history_llr,
                             &level56_shared_invocation,
                             p.LEVEL56_TEMPORAL_LOOKAHEAD_ENABLE
                                 ? &level56_temporal_state
                                 : nullptr,
                             p.LEVEL56_BUFFERED_FIFO_ENABLE
                                 ? &level56_buffered_state
                                 : nullptr);

    win_start += POP_PUSH_ROWS;
  }

  if (p.LEVEL56_BUFFERED_FIFO_ENABLE &&
      p.LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END) {
    drain_level56_buffered_fifo_at_frame_end(
        work_llr, channel_llr, p, stats_ptr, normalize_extrinsic, tx_llr_ref,
        core_fn, &last_tile_history_llr, &level56_shared_invocation,
        &level56_buffered_state);
  }

  if (tile_stats) {
    *tile_stats = std::move(local_tile_stats);
  }

  matrix::Matrix<LLR> out(RROWS, N);
  for (size_t r = 0; r < RROWS; ++r) {
    for (size_t c = 0; c < N; ++c) {
      // 最终输出不是单纯的外信息，而是“原始信道项 + 窗口累积的历史项”。
      const float sum = qfloat::llr_to_float(channel_llr[r][c]) + last_tile_history_llr[r][c];
      out[r][c] = qfloat::llr_from_float<LLR>(sum);
    }
  }

  // 可选：保存窗口累积后的 work_llr（解码前）
  if (p.DUMP_WORK_LLR) {
    matrix::Matrix<float> work_float(RROWS, N);
    for (size_t r = 0; r < RROWS; ++r)
      for (size_t c = 0; c < N; ++c)
        work_float[r][c] = qfloat::llr_to_float(work_llr[r][c]);

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

} // namespace detail
} // namespace newcode
