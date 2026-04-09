#pragma once

#include "newcode/ofec_decoder.hpp"
#include "newcode/params.hpp"
#include "newcode/rx/ofec/chase/chase256.hpp" // 保留 Chase 头；本文档内有三参前向声明
#include "newcode/ofec_decoder_hard.hpp"
#include "newcode/rx/ofec/chase/decoder_core.hpp"
#include "newcode/common/qfloat/llr_utils.hpp"
#include "newcode/decoder_api.hpp"
#include "newcode/ofec/common/lin_matrix_adapters.hpp"
#include "newcode/ofec/earlystop/tile_early_stop_stats.hpp"
#include "newcode/ofec/mux/mux_group_budget.hpp"
#include "newcode/ofec/mux/mux_scheme_c_staged.hpp"
#include "newcode/ofec/mux/mux_state_schedule_apply.hpp"
#include "newcode/ofec/mux/mux_state_builder.hpp"

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
                                                 const std::vector<bool>* early_stop_row_flags,
                                                 const std::vector<uint8_t>* mux_state);

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
                                         int siso_active_for_tile,
                                         bool use_hard_decode,
                                         bool normalize_extrinsic,
                                         const matrix::Matrix<float>* tx_llr_ref,
                                         CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn,
                                         matrix::Matrix<float>* last_tile_history_accum,
                                         bool capture_last_tile_history)
{
  // 输入:
  // - tile_in: 当前 tile 的先验/历史输入矩阵。
  // - ch_tile: 当前 tile 对应的原始信道矩阵。
  // - p: 当前 tile 生效的参数。
  // - tile_top_row_global: tile 顶部在整帧中的全局行号。
  // - siso_active_for_tile: 当前 tile 可用的 SISO budget。
  // - use_hard_decode: 是否走硬判回退。
  // - normalize_extrinsic: 是否对外信息归一化。
  // - tx_llr_ref: 可选参考矩阵，仅调试用。
  // - core_fn: Chase core 回调。
  // - last_tile_history_accum: 可选历史累积矩阵。
  // - capture_last_tile_history: 是否把当前 tile 的历史结果写入 last_tile_history_accum。
  // 输出:
  // - TileProcessResult，包含 tile 写回矩阵以及 early-stop 统计结果。
  // 用途:
  // - 这是 tile 级调度入口，串起“准备输入 -> early-stop -> mux -> core 解码 -> 写回”。
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
      return TileProcessResult<LLR>{tile_out, false, 0, 0, 0, 0};
  }

  TilePrepared<LLR> prep = prepare_tile_inputs(tile_in, ch_tile, p,
                                               tile_top_row_global,
                                               SBR,
                                               rows_to_decode,
                                               tx_llr_ref);

  TileEarlyStopResult early_stop_stats;
  if (p.ENABLE_EARLY_STOP) {
    // 先根据输入统计结果判断哪些 decoder row 已经满足 early-stop 条件。
    early_stop_stats = detect_tile_early_stop(prep.lin_matrix, p);
  } else {
    early_stop_stats.row_passed_flags.assign(rows_to_decode, false);
    early_stop_stats.row_details.assign(rows_to_decode, TileEarlyStopRowDetail{});
    early_stop_stats.rows_passed = 0;
    early_stop_stats.rows_total = rows_to_decode;
    early_stop_stats.all_rows_passed = false;
  }
  bool early_stop_triggered = early_stop_stats.all_rows_passed;
  std::vector<uint8_t> mux_state =
      newcode::mux::build_state_from_early_stop(early_stop_stats);
  const auto count_state =
      [&mux_state](newcode::mux::StateTag tag) -> std::size_t {
        return static_cast<std::size_t>(std::count(
            mux_state.begin(),
            mux_state.end(),
            static_cast<uint8_t>(tag)));
      };
  const std::size_t rows_need_siso_before_mux =
      count_state(newcode::mux::StateTag::NeedSiso);
  // 方案一：只有当当前 tile 的 SISO 预算小于待处理 code 数时，MUX 才真正介入。
  // 若预算已经覆盖全部 rows_to_decode，则保持 NeedSiso/EarlyStopped 原状态，
  // 不再执行 grouped budget 或 reconfig 调度。
  const bool mux_needed =
      siso_active_for_tile < static_cast<int>(rows_to_decode);
  if (mux_needed) {
    const bool use_priority_mux = (p.MUX_SCHEDULING_MODE == 1);
    // 当 SISO 预算不足以覆盖所有 row 时，MUX 才会真正参与裁剪/重配置。
    if (p.MUX_ENABLE_RECONFIG) {
      if (use_priority_mux) {
        newcode::mux::apply_siso_budget_grouped_priority(
            mux_state,
            siso_active_for_tile,
            p.MUX_GROUP_G,
            early_stop_stats,
            p.MUX_EARLY_STOP_PRIORITY_RULE);
      }
      const auto active_codes =
          newcode::mux::collect_active_codes_from_state(mux_state);
      const auto free_siso =
          newcode::mux::build_free_siso_list(siso_active_for_tile);
      const auto schedule = newcode::mux::schedule_scheme_c_staged_cpp(
          static_cast<int>(mux_state.size()),
          siso_active_for_tile,
          p.MUX_GROUP_G,
          active_codes,
          free_siso,
          p.MUX_EXTRA_BYPASS_EDGES);
      newcode::mux::apply_schedule_result_to_mux_state(
          mux_state, schedule.final_code_to_siso);
    } else {
      if (use_priority_mux) {
        newcode::mux::apply_siso_budget_grouped_priority(
            mux_state,
            siso_active_for_tile,
            p.MUX_GROUP_G,
            early_stop_stats,
            p.MUX_EARLY_STOP_PRIORITY_RULE);
      } else {
        newcode::mux::apply_siso_budget_grouped(
            mux_state, siso_active_for_tile, p.MUX_GROUP_G);
      }
    }
  }
  const std::size_t rows_unscheduled =
      count_state(newcode::mux::StateTag::Unscheduled);

  auto decoder_res = decode_tile<LLR>(prep,
                                      use_hard_decode,
                                      normalize_extrinsic,
                                      p,
                                      &early_stop_stats.row_passed_flags,
                                      &mux_state,
                                      core_fn);

  // 将 core 输出重新映射回 tile 的原始坐标系，并可选记录历史值。
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
      early_stop_stats.rows_total,
      rows_need_siso_before_mux,
      rows_unscheduled};
}

} // namespace detail
} // namespace newcode
