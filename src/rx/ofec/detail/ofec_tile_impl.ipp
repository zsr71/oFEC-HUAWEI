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
#include "newcode/ofec/earlystop/tile_early_stop_group_bind.hpp"
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
#include <sstream>
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
#include "ofec_tile_dispatch.ipp"
#include "ofec_tile_decode.ipp"
#include "ofec_tile_writeback.ipp"

namespace newcode {
namespace detail {

inline void increment_hybrid_class_count(HybridRowClass row_class,
                                         HybridClassCount* count) {
  switch (row_class) {
    case HybridRowClass::None:
      ++count->class_none_count;
      return;
    case HybridRowClass::BchHardDecoded:
      ++count->class_bch_hard_decoded_count;
      return;
    case HybridRowClass::Clean:
      ++count->class_clean_count;
      return;
    case HybridRowClass::ParityOnly:
      ++count->class_parity_only_count;
      return;
    case HybridRowClass::OneMain:
      ++count->class_one_main_count;
      return;
    case HybridRowClass::OneMainPlusParity:
      ++count->class_one_main_plus_parity_count;
      return;
    case HybridRowClass::TwoMain:
      ++count->class_two_main_count;
      return;
    case HybridRowClass::Suspicious:
      ++count->class_suspicious_count;
      return;
    case HybridRowClass::HardFail:
      ++count->class_hard_fail_count;
      return;
  }
}

template <typename LLR>
TileProcessResult<LLR> run_legacy_soft_tile_process(
    const matrix::Matrix<LLR>& tile_in,
    const TilePrepared<LLR>& prep,
    const TileEarlyStopResult& early_stop_stats,
    const newcode::Params& p,
    size_t tile_top_row_global,
    int siso_active_for_tile,
    bool normalize_extrinsic,
    CoreFn<typename LinMatrixAdapter<LLR>::core_type> core_fn) {
  // 这是“方案三引入前”的旧 soft tile 路径封装。
  // 只在回归校验模式下调用，用来验证：
  // HYBRID_ENABLE=false 时，新 dispatch-plan 框架是否与旧实现逐 tile 完全一致。
  // 注意：这里不走新的 dispatch plan / hard-finish 分流，只保留旧的 soft tile 语义。
  // 也就是说，后面所有逻辑都要尽量贴近原先的“early-stop -> mux -> decode -> writeback”顺序。
  // 这个函数本身的作用不是改行为，而是给新实现提供一条可对照的旧基线。

  // 先根据 early-stop 结果构造旧路径使用的 mux_state：
  // - NeedSiso：该行还需要软译码资源
  // - EarlyStopped：已命中 early-stop，不再抢 SISO
  // - Unscheduled：后续 MUX 裁掉的行，本轮不产出输出
  // 这个状态数组就是旧 soft path 的核心输入之一。
  std::vector<uint8_t> mux_state =
      newcode::mux::build_state_from_early_stop(early_stop_stats);

  // 用一个小闭包统计 mux_state 里某种状态出现了多少次。
  // 这样后面既能拿到“进入 MUX 前需要 SISO 的行数”，也能拿到“最终被裁掉的行数”。
  const auto count_state =
      [&mux_state](newcode::mux::StateTag tag) -> std::size_t {
        return static_cast<std::size_t>(std::count(
            mux_state.begin(),
            mux_state.end(),
            static_cast<uint8_t>(tag)));
      };

  // 记录进入 MUX 前，真正还在争抢 SISO 的行数。
  // 这条统计是为了和新 dispatch-plan 路径做一致性比较。
  const std::size_t rows_need_siso_before_mux =
      count_state(newcode::mux::StateTag::NeedSiso);

  // 是否真的需要走 MUX 裁剪。
  // 如果 SISO 资源足够，就不必动 mux_state，旧实现也是这样。
  const bool mux_needed =
      siso_active_for_tile < static_cast<int>(prep.lin_matrix.rows());
  if (mux_needed) {
    // 旧路径同样支持两种预算分配方式：
    // - 普通 group budget
    // - 带 early-stop 优先级的 group budget
    // 这里保留与新路径一致的参数选择，只是不引入新的 dispatch plan 结构。
    const bool use_priority_mux = (p.MUX_SCHEDULING_MODE == 1);
    if (p.MUX_ENABLE_RECONFIG) {
      // reconfig 模式下，旧路径也会走 scheme-C staged 调度。
      // 这里仍然是直接在 mux_state 上原地改写，不做额外的 plan 封装。
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
      // staged 调度器会根据当前活跃 code 和空闲 SISO 槽位，重新生成 code->SISO 映射。
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

  // 统计 MUX 结束后仍然被标成 Unscheduled 的行数。
  // 这是旧路径的一个关键结果，也要和新路径的统计对齐。
  const std::size_t rows_unscheduled =
      count_state(newcode::mux::StateTag::Unscheduled);

  // 旧路径直接调用 tile decoder：
  // - use_hard_decode=false：走 soft/Chase 路径
  // - early_stop_row_flags：把 row early-stop 结果传下去
  // - mux_state：告诉 core 哪些行真正参与 soft decode，哪些行要跳过
  // 这一步是旧 soft tile 语义的核心执行点。
  auto decoder_res = decode_tile<LLR>(prep,
                                      /*use_hard_decode=*/false,
                                      normalize_extrinsic,
                                      p,
                                      &early_stop_stats.row_passed_flags,
                                      &mux_state,
                                      core_fn);

  // 旧路径的结果写回也保持原样：
  // 先在 tile_out 上做回填，再把最终 tile 输出交给 TileProcessResult。
  matrix::Matrix<LLR> tile_out = tile_in;
  writeback_tile(prep,
                 decoder_res,
                 p,
                 tile_top_row_global,
                 /*capture_last_tile_history=*/false,
                 &tile_out,
                 /*last_tile_history_accum=*/nullptr);

  // 返回值也沿用旧路径的统计口径：
  // - early-stop 是否触发
  // - early-stop 行数
  // - total 行数
  // - hard_finish 固定为 0（旧路径没有 hybrid prepass）
  // - mux 前 need_siso 数
  // - MUX 后 unscheduled 数
  return TileProcessResult<LLR>{
      std::move(tile_out),
      early_stop_stats.all_rows_passed,
      early_stop_stats.rows_passed,
      early_stop_stats.rows_total,
      /*rows_hard_finish=*/0,
      rows_need_siso_before_mux,
      rows_unscheduled};
}

template <typename LLR>
void validate_hybrid_disabled_matches_legacy(
    const TileProcessResult<LLR>& actual,
    const TileProcessResult<LLR>& legacy,
    size_t tile_top_row_global) {
  // 这个函数只在 `HYBRID_ENABLE=false` 的回归模式下使用。
  // 它的职责不是修复差异，而是确认“新 dispatch-plan 路径”与“旧 soft 路径”逐 tile 行为一致。
  // 只要发现差异，就立刻抛错，让问题停在第一次出现的 tile 上，方便定位。

  // 下面的 fail 闭包负责统一拼接错误信息：
  // - tile_top_row_global：告诉你是帧里的哪一个 tile
  // - what：告诉你具体是哪一种字段或哪一个坐标出现了不一致
  auto fail = [&](const std::string& what) -> void {
    std::ostringstream oss;
    oss << "HYBRID_ENABLE=false regression mismatch at tile_top_row_global="
        << tile_top_row_global << ": " << what;
    throw std::runtime_error(oss.str());
  };

  // 先比较的是“是否触发 early-stop”这个最外层状态。
  // 这一步如果不一致，说明新旧路径在最上游分流就已经偏了。
  if (actual.early_stop_triggered != legacy.early_stop_triggered) {
    fail("early_stop_triggered differs");
  }
  // 再比较 early-stop 相关统计量，确保两条路径对行级 early-stop 的判定口径一致。
  if (actual.rows_early_stop != legacy.rows_early_stop) {
    fail("rows_early_stop differs");
  }
  if (actual.rows_total != legacy.rows_total) {
    fail("rows_total differs");
  }
  // `rows_hard_finish` 在 hybrid 关闭时应当恒等于 0；
  // 这里依然比较，是为了防止新路径把统计口径偷偷改掉。
  if (actual.rows_hard_finish != legacy.rows_hard_finish) {
    fail("rows_hard_finish differs");
  }
  // 这项比较的是 MUX 前、真正需要 SISO 的行数。
  // 它能反映新旧路径在“候选集合”构建上的一致性。
  if (actual.rows_need_siso_before_mux != legacy.rows_need_siso_before_mux) {
    fail("rows_need_siso_before_mux differs");
  }
  // 这项比较的是 MUX 后被裁掉的行数。
  // 如果这里不同，通常说明预算分配或状态回写已经偏离旧实现。
  if (actual.rows_unscheduled != legacy.rows_unscheduled) {
    fail("rows_unscheduled differs");
  }
  // 先比较输出矩阵形状，防止后面的逐元素比较在维度不一致时误判。
  if (actual.tile_out.rows() != legacy.tile_out.rows() ||
      actual.tile_out.cols() != legacy.tile_out.cols()) {
    fail("tile_out shape differs");
  }

  // 最后逐元素比较 tile_out。
  // 这里用一个很小的浮点容差，是为了避免仅由量化/打印误差引入的假阳性。
  constexpr float kTol = 1e-6f;
  for (size_t r = 0; r < actual.tile_out.rows(); ++r) {
    for (size_t c = 0; c < actual.tile_out.cols(); ++c) {
      const float a = qfloat::llr_to_float(actual.tile_out[r][c]);
      const float b = qfloat::llr_to_float(legacy.tile_out[r][c]);
      if (std::fabs(a - b) > kTol) {
        std::ostringstream oss;
        oss << "tile_out differs at (" << r << "," << c
            << "), actual=" << a
            << ", legacy=" << b;
        fail(oss.str());
      }
    }
  }
}

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
    // 条件1可额外按固定 group 绑定：只有整组 row 都通过时，这组才整体 early-stop。
    if (p.EARLY_STOP_CONDITION_MODE == 1 &&
        p.EARLY_STOP_BIND_GROUP_SIZE > 1) {
      early_stop_stats = apply_group_bound_early_stop(
          early_stop_stats, p.EARLY_STOP_BIND_GROUP_SIZE);
    }
  } else {
    early_stop_stats.row_passed_flags.assign(rows_to_decode, false);
    early_stop_stats.row_details.assign(rows_to_decode, TileEarlyStopRowDetail{});
    early_stop_stats.rows_passed = 0;
    early_stop_stats.rows_total = rows_to_decode;
    early_stop_stats.all_rows_passed = false;
  }
  bool early_stop_triggered = early_stop_stats.all_rows_passed;
  std::size_t rows_need_siso_before_mux = 0;
  std::size_t rows_hard_finish = 0;
  std::size_t rows_unscheduled = 0;
  HybridClassCount hybrid_class_count{};
  const bool verify_hybrid_disabled =
      !use_hard_decode &&
      p.HYBRID_VERIFY_DISABLED_MATCH_LEGACY;
  // 这两个限制是为了保证比较语义纯粹：
  // - hybrid 必须真的关闭
  // - normalize 语义必须保持旧行为
  if (verify_hybrid_disabled && p.HYBRID_ENABLE) {
    throw std::invalid_argument(
        "HYBRID_VERIFY_DISABLED_MATCH_LEGACY requires HYBRID_ENABLE=false.");
  }
  if (verify_hybrid_disabled && p.HYBRID_NORMALIZE_SOFT_ONLY) {
    throw std::invalid_argument(
        "HYBRID_VERIFY_DISABLED_MATCH_LEGACY requires HYBRID_NORMALIZE_SOFT_ONLY=false.");
  }

  chase::DecoderCoreResult<typename LinMatrixAdapter<LLR>::core_type> decoder_res;
  if (use_hard_decode) {
    std::vector<uint8_t> mux_state =
        newcode::mux::build_state_from_early_stop(early_stop_stats);
    const auto count_state =
        [&mux_state](newcode::mux::StateTag tag) -> std::size_t {
          return static_cast<std::size_t>(std::count(
              mux_state.begin(),
              mux_state.end(),
              static_cast<uint8_t>(tag)));
        };
    rows_need_siso_before_mux =
        count_state(newcode::mux::StateTag::NeedSiso);
    const bool mux_needed =
        siso_active_for_tile < static_cast<int>(rows_to_decode);
    if (mux_needed) {
      const bool use_priority_mux = (p.MUX_SCHEDULING_MODE == 1);
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
    rows_unscheduled = count_state(newcode::mux::StateTag::Unscheduled);

    decoder_res = decode_tile<LLR>(prep,
                                   use_hard_decode,
                                   normalize_extrinsic,
                                   p,
                                   &early_stop_stats.row_passed_flags,
                                   &mux_state,
                                   core_fn);
  } else {
    // soft tile 走方案三新框架：
    // 1. 建 dispatch plan
    // 2. 只让 soft candidate 参加 MUX
    // 3. 按计划执行并合并结果
    auto dispatch_plan = build_tile_dispatch_plan(
        prep, early_stop_stats, siso_active_for_tile, p);
    rows_need_siso_before_mux = count_soft_candidates_before_mux(dispatch_plan);
    rows_hard_finish = dispatch_plan.rows_hard_finish;
    hybrid_class_count.invocation = static_cast<std::size_t>(
        std::max(0, p.debug_trace.chase_invocation));
    hybrid_class_count.tile_index = static_cast<std::size_t>(
        std::max(0, p.debug_trace.chase_tile_index));
    hybrid_class_count.rows_seen_by_hybrid = dispatch_plan.rows_seen_by_hybrid;
    hybrid_class_count.deferred_candidate_count =
        dispatch_plan.deferred_candidate_count;
    hybrid_class_count.deferred_priority_0_count =
        dispatch_plan.deferred_priority_0_count;
    hybrid_class_count.deferred_priority_1_count =
        dispatch_plan.deferred_priority_1_count;
    hybrid_class_count.deferred_priority_2_count =
        dispatch_plan.deferred_priority_2_count;
    hybrid_class_count.deferred_priority_3_count =
        dispatch_plan.deferred_priority_3_count;
    hybrid_class_count.deferred_reclaimed_to_hard_finish_count =
        dispatch_plan.deferred_reclaimed_to_hard_finish_count;
    for (std::size_t row = 0; row < dispatch_plan.hybrid_classes.size(); ++row) {
      if (row < dispatch_plan.rows.size() &&
          dispatch_plan.rows[row].early_stop_hit) {
        continue;
      }
      increment_hybrid_class_count(dispatch_plan.hybrid_classes[row],
                                   &hybrid_class_count);
    }
    run_mux_on_soft_candidates(&dispatch_plan,
                               siso_active_for_tile,
                               p,
                               early_stop_stats);
    rows_unscheduled = dispatch_plan.rows_soft_unscheduled;
    decoder_res = decode_tile_with_plan<LLR>(prep,
                                             dispatch_plan,
                                             normalize_extrinsic,
                                             core_fn);
  }

  // 将 core 输出重新映射回 tile 的原始坐标系，并可选记录历史值。
  writeback_tile(prep,
                 decoder_res,
                 p,
                 tile_top_row_global,
                 capture_last_tile_history,
                 &tile_out,
                 last_tile_history_accum);

  TileProcessResult<LLR> tile_result{
      std::move(tile_out),
      early_stop_triggered,
      early_stop_stats.rows_passed,
      early_stop_stats.rows_total,
      rows_hard_finish,
      rows_need_siso_before_mux,
      rows_unscheduled,
      std::move(hybrid_class_count)};

  if (verify_hybrid_disabled) {
    // 调试校验模式下，同一块 tile 再跑一遍旧 soft 路径做逐项比对。
    const auto legacy_result = run_legacy_soft_tile_process(
        tile_in,
        prep,
        early_stop_stats,
        p,
        tile_top_row_global,
        siso_active_for_tile,
        normalize_extrinsic,
        core_fn);
    validate_hybrid_disabled_matches_legacy(
        tile_result, legacy_result, tile_top_row_global);
  }

  return tile_result;
}

} // namespace detail
} // namespace newcode
