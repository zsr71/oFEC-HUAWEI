#pragma once

#include "newcode/ofec/hybrid/hybrid_classifier.hpp"

#include <algorithm>
#include <array>
#include <vector>

namespace newcode {
namespace detail {

// 方案三里把“资源竞争语义”和“执行动作语义”拆开后，
// 每一行最终只会落到这四种执行标签之一。
enum class RowDispatchTag : uint8_t {
  SoftDecode = 0,
  EarlyStopAction,
  HardFinish,
  Unscheduled
};

// 逐行调度主表：
// - tag: 这行最终走哪条执行路径
// - hybrid_class: hybrid prepass 给出的原因标签
// - early_stop_hit: 这行是否在 row early-stop 判定中命中
// - scheduled_for_soft: 在参与 MUX 后，这行是否真的抢到了 SISO
struct RowDispatchEntry {
  RowDispatchTag tag = RowDispatchTag::SoftDecode;
  HybridRowClass hybrid_class = HybridRowClass::None;
  bool early_stop_hit = false;
  bool scheduled_for_soft = false;
};

template <typename CoreLLR>
struct TileDispatchPlan {
  // 完整行域的调度主表。rows[i] 对应第 i 个 decoder row。
  std::vector<RowDispatchEntry> rows;

  // 这三张表描述 soft path 的筛选过程：
  // candidate -> scheduled / unscheduled。
  std::vector<int> soft_candidate_rows;
  std::vector<int> soft_scheduled_rows;
  std::vector<int> soft_unscheduled_rows;

  // 前置硬纠成功时，直接把该行的 lout 缓存在这里；
  // decode 阶段无需再跑 Chase，只需把结果回填。
  matrix::Matrix<float> hard_finish_lout;
  std::vector<bool> hard_finish_valid;

  // 调试/统计专用的分类结果镜像，便于窗口级汇总和后续扩展。
  std::vector<HybridRowClass> hybrid_classes;

  // 各类行数摘要，避免后面重复扫描 rows。
  std::size_t rows_early_stop = 0;
  std::size_t rows_hard_finish = 0;
  std::size_t rows_soft_candidate = 0;
  std::size_t rows_soft_scheduled = 0;
  std::size_t rows_soft_unscheduled = 0;
};

template <typename LLR>
inline void load_row_vectors(const TilePrepared<LLR>& prep,
                             std::size_t row,
                             std::array<typename TilePrepared<LLR>::CoreLLR,
                                        newcode::Params::BCH_N>* lin_vec,
                             std::array<typename TilePrepared<LLR>::CoreLLR,
                                        newcode::Params::BCH_N>* lch_vec) {
  // 从矩阵视图里拷出一整行，便于复用现有的行级硬纠接口。
  for (std::size_t col = 0; col < static_cast<std::size_t>(newcode::Params::BCH_N);
       ++col) {
    (*lin_vec)[col] = prep.lin_matrix[row][col];
    (*lch_vec)[col] = prep.lch_matrix[row][col];
  }
}

template <typename LLR>
void rebuild_soft_candidate_rows(
    TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>* plan) {
  // 有些行在 hybrid prepass 后会从 SoftDecode 改成 HardFinish，
  // 所以需要根据 rows 主表重新生成一次 soft 候选名单。
  plan->soft_candidate_rows.clear();  // 先清空旧的候选行列表，避免混入上一轮/上一阶段的残留结果。

  // 遍历完整行域（decoder row 维度），逐行判断“是否仍需要参与 soft decode”。
  for (std::size_t row = 0; row < plan->rows.size(); ++row) {
    // 只有 tag==SoftDecode 的行才是 soft 候选：
    // - EarlyStopAction：已早停，不进 Chase
    // - HardFinish：已硬纠完成，不进 Chase
    // - Unscheduled：被 MUX 裁掉，本轮不进 Chase
    if (plan->rows[row].tag == RowDispatchTag::SoftDecode) {
      // 记录 full row index；后续 MUX 和 soft-only mux_state 会直接用这些原始行号。
      plan->soft_candidate_rows.push_back(static_cast<int>(row));
    }
  }

  // 同步候选计数，供 tile/window 统计与日志输出使用（避免后面再扫一次 vector）。
  plan->rows_soft_candidate = plan->soft_candidate_rows.size();
}

inline bool hybrid_class_is_siso_backfill_candidate(
    HybridRowClass row_class,
    newcode::HybridSisoBackfillMode mode) {
  switch (mode) {
    case newcode::HybridSisoBackfillMode::TwoErrorOnly:
      return row_class == HybridRowClass::TwoMain;
    case newcode::HybridSisoBackfillMode::Disabled:
    default:
      return false;
  }
}

template <typename CoreLLR>
inline void apply_hard_finish_to_plan(
    std::size_t row,
    HybridRowClass hard_class,
    const std::array<float, newcode::Params::BCH_N>& y2,
    TileDispatchPlan<CoreLLR>* plan) {
  plan->rows[row].tag = RowDispatchTag::HardFinish;
  plan->rows[row].hybrid_class = hard_class;
  plan->hybrid_classes[row] = hard_class;
  plan->hard_finish_valid[row] = true;
  for (std::size_t col = 0; col < static_cast<std::size_t>(newcode::Params::BCH_N);
       ++col) {
    plan->hard_finish_lout[row][col] = y2[col];
  }
  ++plan->rows_hard_finish;
}

template <typename LLR>
void run_hybrid_prepass(
    const TilePrepared<LLR>& prep,
    const std::vector<bool>& early_stop_row_flags,
    int siso_active_for_tile,
    const newcode::Params& p,
    TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR>* plan) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  if (!p.HYBRID_ENABLE) {
    // 关闭 hybrid 时，plan 仍然统一走方案三框架，
    // 但此时不做任何前置分流，保持所有非 early-stop 行继续走 soft path。
    rebuild_soft_candidate_rows<LLR>(plan);
    return;
  }

  // hybrid prepass 的目标：
  // - 对“尚未 early-stop 且仍处于 SoftDecode”的行，先尝试一把硬纠（BCH hard decode）
  // - 硬纠成功：通常会变成 HardFinish，本轮直接产出输出，不再进入 Chase
  // - 当开启 SISO backfill 时，部分 hard-finish 候选可以被保留给后续 soft path，以尽量用满 SISO
  // - 硬纠失败：这行保留 SoftDecode，后续仍参与 MUX/SISO 竞争与 Chase
  const auto classifier_mode = effective_hybrid_classifier_mode(p);
  const bool enable_siso_backfill =
      classifier_mode != newcode::HybridClassifierMode::LegacyHardDecode &&
      p.HYBRID_SISO_BACKFILL_MODE != newcode::HybridSisoBackfillMode::Disabled;
  struct DeferredHardFinishCandidate {
    std::size_t row = 0;
    HybridRowClass hard_class = HybridRowClass::None;
    std::array<float, newcode::Params::BCH_N> y2{};
  };
  std::vector<DeferredHardFinishCandidate> deferred_candidates;
  const std::size_t rows = prep.lin_matrix.rows();
  for (std::size_t row = 0; row < rows; ++row) {
    if (row < early_stop_row_flags.size() && early_stop_row_flags[row]) {
      // 已经命中 early-stop 的行不会再进 hybrid prepass。
      continue;
    }
    if (plan->rows[row].tag != RowDispatchTag::SoftDecode) {
      // 只对 SoftDecode 行做 prepass。
      // EarlyStopAction / HardFinish / Unscheduled 都已经确定走向，不应再次改写。
      continue;
    }

    // 从 tile 级矩阵视图中拷出一整行，便于调用现有行级硬纠接口。
    std::array<CoreLLR, newcode::Params::BCH_N> lin_vec{};
    std::array<CoreLLR, newcode::Params::BCH_N> lch_vec{};
    std::array<float, newcode::Params::BCH_N> y2{};
    load_row_vectors(prep, row, &lin_vec, &lch_vec);

    // hybrid prepass 引擎支持四种模式：
    // - LegacyHardDecode：保持原 E1，直接调用完整 BCH 硬译码
    // - RepoFastClassifier：当前仓库的 S0/S1/S3 快速分类
    // - FriendS1S3Classifier：朋友版 S1/S3 + Trace(mu) 思路
    // - FriendS1S3WithS0Classifier：朋友版 + S0 约束的 2 错入口
    HybridRowClass hard_class = HybridRowClass::HardFail;
    bool hard_ok = false;
    if (classifier_mode == newcode::HybridClassifierMode::LegacyHardDecode) {
      // 兼容旧路径：不经过 fast classifier，直接调用原完整 BCH 硬译码。
      hard_ok = newcode::perform_hard_decode<CoreLLR>(lin_vec, lch_vec, y2, p);
      hard_class = hard_ok ? HybridRowClass::BchHardDecoded
                           : HybridRowClass::HardFail;
    } else {
      // 新路径：先进入独立的分类器模块，由它决定：
      // - 是否能直接 hard-finish
      // - 如果能，属于哪种 hard-finish 分类
      // - 如果不能，失败原因是什么
      hard_ok = run_selected_hybrid_classifier_hard_finish(
          lin_vec, &y2, p, &hard_class);
    }
    if (!hard_ok) {
      // 标记“未进入 hard-finish”的原因标签：
      // - HardFail: 明确不满足当前 hard path 条件
      // - Suspicious: 预留给“分类上可疑，但不能直接产出 hard-finish”的情况
      //
      // 注意：这里不会把 tag 改成 Unscheduled。
      // 也就是说，prepass 失败并不等于这一行被丢弃，
      // 它仍然保留 SoftDecode 身份，后面继续参加 MUX / SISO 竞争。
      plan->rows[row].hybrid_class = hard_class;
      plan->hybrid_classes[row] = hard_class;
      continue;
    }

    if (enable_siso_backfill &&
        hybrid_class_is_siso_backfill_candidate(
            hard_class, p.HYBRID_SISO_BACKFILL_MODE)) {
      // 这类行在分类上已经可以 hard-finish，但当前策略允许先保留为 soft 候选，
      // 等整块 tile 的分类结果都出来后，再结合 SISO 预算决定要不要真正收掉。
      plan->rows[row].hybrid_class = hard_class;
      plan->hybrid_classes[row] = hard_class;
      deferred_candidates.push_back(
          DeferredHardFinishCandidate{row, hard_class, y2});
      continue;
    }

    // 硬纠成功：将该行从 SoftDecode 改为 HardFinish，并缓存输出。
    // 从这一刻起，这行不再进入 Chase；decode 阶段只需要把缓存的 y2 直接回填。
    apply_hard_finish_to_plan(row, hard_class, y2, plan);
  }

  if (enable_siso_backfill && !deferred_candidates.empty()) {
    // 当前策略先把 deferred 候选都留在 SoftDecode 集合里。
    // 如果此时剩余 soft 行数仍然多于 SISO 预算，说明即便收掉一部分 deferred，
    // 也不会导致 SISO 吃不满；这种情况下再把多出来的 deferred 收回 HardFinish。
    const std::size_t soft_rows_before_deferred_accept = static_cast<std::size_t>(
        std::count_if(plan->rows.begin(),
                      plan->rows.end(),
                      [](const RowDispatchEntry& entry) {
                        return entry.tag == RowDispatchTag::SoftDecode;
                      }));
    const std::size_t siso_budget =
        static_cast<std::size_t>(std::max(siso_active_for_tile, 0));
    const std::size_t max_deferred_accept =
        (soft_rows_before_deferred_accept > siso_budget)
            ? (soft_rows_before_deferred_accept - siso_budget)
            : 0u;
    const std::size_t deferred_accept_count =
        std::min(max_deferred_accept, deferred_candidates.size());
    for (std::size_t i = 0; i < deferred_accept_count; ++i) {
      apply_hard_finish_to_plan(deferred_candidates[i].row,
                                deferred_candidates[i].hard_class,
                                deferred_candidates[i].y2,
                                plan);
    }
  }

  // prepass 会把一部分 SoftDecode 行改成 HardFinish，
  // 所以需要在结束时重建一次 soft candidate 列表，保证后续 MUX 输入正确。
  rebuild_soft_candidate_rows<LLR>(plan);
}

template <typename LLR>
TileDispatchPlan<typename TilePrepared<LLR>::CoreLLR> build_tile_dispatch_plan(
    const TilePrepared<LLR>& prep,
    const TileEarlyStopResult& early_stop_stats,
    int siso_active_for_tile,
    const newcode::Params& p) {
  using CoreLLR = typename TilePrepared<LLR>::CoreLLR;
  TileDispatchPlan<CoreLLR> plan;
  const std::size_t rows = prep.lin_matrix.rows();
  const std::size_t cols = prep.lin_matrix.cols();

  // 第 1 步：按完整 decoder row 域初始化调度主表。
  // 这里先把所有 row / 缓存表都分配好，后续只做逐行改写，不再反复扩容。
  plan.rows.assign(rows, RowDispatchEntry{});
  // hybrid_classes 是 debug / 统计镜像表，先全部置为 None，后面再按实际分类写入。
  plan.hybrid_classes.assign(rows, HybridRowClass::None);
  // hard_finish_lout 用来缓存 hybrid prepass 的硬纠结果，按完整行域分配。
  plan.hard_finish_lout = matrix::Matrix<float>(rows, cols);
  // hard_finish_valid 标记每个 row 的缓存是否真的有效，初始全部 false。
  plan.hard_finish_valid.assign(rows, false);

  // 第 2 步：把 early-stop 判定结果写进主表。
  // 这一步只做“是否早停”的分流，不做 soft/hard 的资源竞争判断。
  // 第一步先把完整行域切成：
  // - early-stop 行
  // - 其余待定行（先暂标成 SoftDecode）
  for (std::size_t row = 0; row < rows; ++row) {
    // 安全读取每一行的 early-stop 命中结果；越界时视为 false。
    const bool early_stop_hit =
        row < early_stop_stats.row_passed_flags.size() &&
        early_stop_stats.row_passed_flags[row];
    // 先记录这个布尔值，便于后续调试判断该行是否曾命中过 early-stop。
    plan.rows[row].early_stop_hit = early_stop_hit;
    if (early_stop_hit) {
      // 命中 early-stop 的行直接走 early-stop action，不再参与后续 MUX / soft decode。
      plan.rows[row].tag = RowDispatchTag::EarlyStopAction;
      ++plan.rows_early_stop;
    } else {
      // 未命中的行先作为 SoftDecode 候选，等待 hybrid prepass 或 MUX 再细分。
      plan.rows[row].tag = RowDispatchTag::SoftDecode;
      plan.soft_candidate_rows.push_back(static_cast<int>(row));
    }
  }
  // 这里先保留一份“hybrid prepass 之前”的 soft 候选数，供统计和日志使用。
  plan.rows_soft_candidate = plan.soft_candidate_rows.size();

  // 第 3 步：运行 hybrid prepass。
  // 关闭 HYBRID_ENABLE 时，这一步不会改变行的调度语义，只会重建 soft 候选列表；
  // 打开时，部分 SoftDecode 行会被升级成 HardFinish。
  // 第二步再让 hybrid prepass 把一部分 SoftDecode 行拿走变成 HardFinish。
  run_hybrid_prepass(
      prep, early_stop_stats.row_passed_flags, siso_active_for_tile, p, &plan);
  // 返回完整 tile 的调度计划：
  // - 早停行：EarlyStopAction
  // - hybrid 成功行：HardFinish
  // - 剩余行：SoftDecode，等待后续 MUX / soft core 处理
  return plan;
}

template <typename CoreLLR>
std::size_t count_soft_candidates_before_mux(
    const TileDispatchPlan<CoreLLR>& plan) {
  // 这里的计数语义是：
  // “经过 early-stop + hybrid prepass 之后，还剩多少行要去竞争 SISO”。
  return plan.soft_candidate_rows.size();
}

template <typename CoreLLR>
void run_mux_on_soft_candidates(TileDispatchPlan<CoreLLR>* plan,
                                int siso_active_for_tile,
                                const newcode::Params& p,
                                const TileEarlyStopResult& early_stop_stats) {
  // 先把上一次调度留下来的结果清掉，避免复用旧 plan 时把历史状态带进来。
  plan->soft_scheduled_rows.clear();
  plan->soft_unscheduled_rows.clear();

  // 为了最大程度复用现有 MUX 实现，这里构造一个“伪 mux_state”：
  // - 只有 SoftDecode 行标成 NeedSiso
  // - 其他行一律视作已退出竞争（EarlyStopped 占位）
  std::vector<uint8_t> mux_candidate_state(
      plan->rows.size(),
      static_cast<uint8_t>(newcode::mux::StateTag::EarlyStopped));
  for (std::size_t row = 0; row < plan->rows.size(); ++row) {
    // 只要这行在 dispatch plan 里仍是 SoftDecode，就让它进入 SISO 竞争池。
    // 其他 tag 的行都不应该再参与 MUX，因此保持 EarlyStopped 占位即可。
    if (plan->rows[row].tag == RowDispatchTag::SoftDecode) {
      mux_candidate_state[row] =
          static_cast<uint8_t>(newcode::mux::StateTag::NeedSiso);
    }
  }

  // 统计当前还有多少 soft 候选行真正要去争抢 SISO。
  // 这个数就是后续判断“需不需要裁剪”的依据。
  const auto need_count = static_cast<int>(std::count(
      mux_candidate_state.begin(),
      mux_candidate_state.end(),
      static_cast<uint8_t>(newcode::mux::StateTag::NeedSiso)));

  // 只有 soft 候选数超过本 tile 的 SISO 预算时，才真正触发裁剪。
  if (need_count > 0 && siso_active_for_tile < need_count) {
    // MUX 有两层可选策略：
    // - MUX_SCHEDULING_MODE==0: 按常规预算分配
    // - MUX_SCHEDULING_MODE==1: 先按 early-stop 结果做优先级调度，再分配预算
    const bool use_priority_mux = (p.MUX_SCHEDULING_MODE == 1);
    if (p.MUX_ENABLE_RECONFIG) {
      // reconfig 版本会走 scheme-C staged 流程：
      // 先可选做优先级打分，再把活跃 code 和 free SISO 列表喂给 staged 调度器。
      if (use_priority_mux) {
        newcode::mux::apply_siso_budget_grouped_priority(
            mux_candidate_state,
            siso_active_for_tile,
            p.MUX_GROUP_G,
            early_stop_stats,
            p.MUX_EARLY_STOP_PRIORITY_RULE);
      }
      // 从当前 state 里收集哪些 code 还活着，作为 staged 调度的输入。
      const auto active_codes =
          newcode::mux::collect_active_codes_from_state(mux_candidate_state);
      // free_siso 表示本 tile 里可用的 SISO 资源槽位列表。
      const auto free_siso =
          newcode::mux::build_free_siso_list(siso_active_for_tile);
      // staged 调度器会根据当前活跃 code、空闲 SISO 和额外 bypass 边，
      // 重新构造 code -> siso 的映射。
      const auto schedule = newcode::mux::schedule_scheme_c_staged_cpp(
          static_cast<int>(mux_candidate_state.size()),
          siso_active_for_tile,
          p.MUX_GROUP_G,
          active_codes,
          free_siso,
          p.MUX_EXTRA_BYPASS_EDGES);
      // 把 staged 调度结果写回 state：能拿到 SISO 的行保留 NeedSiso，拿不到的行会被清掉。
      newcode::mux::apply_schedule_result_to_mux_state(
          mux_candidate_state, schedule.final_code_to_siso);
    } else {
      // 非 reconfig 路径直接用现有预算分配器。
      if (use_priority_mux) {
        // 优先级模式会把 early-stop 命中情况纳入排序依据。
        newcode::mux::apply_siso_budget_grouped_priority(
            mux_candidate_state,
            siso_active_for_tile,
            p.MUX_GROUP_G,
            early_stop_stats,
            p.MUX_EARLY_STOP_PRIORITY_RULE);
      } else {
        // 普通模式按 group 预算直接裁剪，不引入额外优先级规则。
        newcode::mux::apply_siso_budget_grouped(
            mux_candidate_state, siso_active_for_tile, p.MUX_GROUP_G);
      }
    }
  }

  // 把 MUX 的结果重新写回 dispatch plan 主表。
  for (std::size_t row = 0; row < plan->rows.size(); ++row) {
    // 只有 soft 候选行需要根据 MUX 结果改写；
    // 其他 tag 的行在前面已经确定，不应该被这里再次修改。
    if (plan->rows[row].tag != RowDispatchTag::SoftDecode) {
      continue;
    }
    // 仍然保留 NeedSiso 的行，说明它成功抢到了预算，后面会进入 soft decode。
    if (mux_candidate_state[row] ==
        static_cast<uint8_t>(newcode::mux::StateTag::NeedSiso)) {
      plan->rows[row].scheduled_for_soft = true;
      plan->soft_scheduled_rows.push_back(static_cast<int>(row));
    } else {
      // 没抢到预算的 soft 候选，本轮直接标成 Unscheduled。
      // 这样 decode 阶段构造 soft-only mux_state 时就不会让它跑 Chase。
      plan->rows[row].tag = RowDispatchTag::Unscheduled;
      plan->rows[row].scheduled_for_soft = false;
      plan->soft_unscheduled_rows.push_back(static_cast<int>(row));
    }
  }

  // 同步统计值，供 tile/window/pipeline 汇总输出。
  plan->rows_soft_scheduled = plan->soft_scheduled_rows.size();
  plan->rows_soft_unscheduled = plan->soft_unscheduled_rows.size();
}

} // namespace detail
} // namespace newcode
