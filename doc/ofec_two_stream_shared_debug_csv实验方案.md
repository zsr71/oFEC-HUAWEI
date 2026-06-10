# `two-stream shared` 调试 CSV 实验方案

本文说明：

- 为了给 `two-stream shared` 增加更完整的内部诊断能力
- 需要改哪些代码
- 每个文件准备怎么改
- 这些改动分几步推进
- 每一步改完后如何验证

这份方案对应的目标，不是修改 shared 解码算法本身，而是：

> 在不改变当前 shared 解码语义的前提下，把内部已经存在或容易获得的调度/分类/结果信息稳定导出来。

也就是说，第一版工作的重点是：

- 增加可观测性
- 增加 probe / CSV 导出能力
- 不主动改 shared 调度决策
- 不主动改 Chase 算法逻辑

---

## 1. 目标定义

这次实验的直接目标是落地三份调试 CSV：

1. `per_invocation_shared_tile_samples.csv`
2. `per_invocation_shared_hybrid_class_counts.csv`
3. `per_invocation_shared_row_map.csv`

它们的口径说明已经在：

- [ofec_two_stream_shared_debug_csv设计说明.md](/home/zsr71/projects/newcode_two_stream_shared_chase/doc/ofec_two_stream_shared_debug_csv设计说明.md:1)

本文不再重复解释列含义，而是回答：

- 这些数据从哪里来
- 哪些代码已经有现成基础
- 哪些地方还要补结构
- 最小侵入的实现路径是什么

---

## 2. 当前代码基础

当前仓库里已经有三类现成基础，可以直接复用。

### 2.1 shared runner 已经能拿到 tile 级摘要

在
[src/two_stream_shared/two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:430)
的 `run_shared_tile(...)` 里，当前已经可以拿到：

- `rows_early_stop`
- `rows_total`
- `rows_hard_finish`
- `rows_need_siso_before_mux`
- `rows_unscheduled`
- `produced_rows`
- `failed_rows`
- `produced_rows_a`
- `produced_rows_b`
- `failed_rows_a`
- `failed_rows_b`

这意味着：

- 第一份 CSV 的大部分列，已经不需要再深入 Chase 核心去新算一遍
- 只需要把当前函数内部已有的统计结果正式挂出来

### 2.2 dispatch plan 已经有 row 级语义

在
[src/rx/ofec/detail/ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:12)
里，`TileDispatchPlan` 和 `RowDispatchEntry` 当前已经保存了：

- `tag`
  - `SoftDecode`
  - `EarlyStopAction`
  - `HardFinish`
  - `Unscheduled`
- `hybrid_class`
- `early_stop_hit`
- `scheduled_for_soft`

这意味着：

- 第三份 CSV 的核心 row 去向信息，本质上已经存在于内存里
- 目前的问题不是“没有信息”，而是“没有导出成稳定结构”

### 2.3 hybrid 分类枚举已经存在

在
[include/newcode/ofec/hybrid/hybrid_classifier.hpp](/home/zsr71/projects/newcode_two_stream_shared_chase/include/newcode/ofec/hybrid/hybrid_classifier.hpp:11)
里，`HybridRowClass` 已经定义了：

- `None`
- `BchHardDecoded`
- `Clean`
- `ParityOnly`
- `OneMain`
- `OneMainPlusParity`
- `TwoMain`
- `Suspicious`
- `HardFail`

这意味着：

- 第二份 CSV 的分类统计，不需要重新定义新的分类系统
- 只需要把当前 row 级分类做计数汇总

---

## 3. 实现总策略

推荐采用下面这条路径：

### 3.1 不直接把 CSV 落盘逻辑塞进 shared runner

shared runner 的职责应保持为：

- 跑两路前端
- 进入 shared merged tile 路径
- 返回主结果和 observability

不建议让它直接承担：

- 文件管理
- 多 seed 调度
- CSV 表头控制
- probe 输出目录组织

因此推荐的职责分层是：

1. shared runner 内部负责“生产结构化调试结果”
2. 新增一个 two-stream probe app 负责“组织任务并落盘 CSV”

### 3.2 第一版优先做“结构化结果导出”，再做 probe app

推荐顺序是：

1. 先把 shared runner 的 observability 结构补齐
2. 确认 one-shot app 能打印或拿到这些调试结果
3. 再新增 `ofec_two_stream_ber_window_probe.cpp`
4. 最后由 probe app 写出三份 CSV

这样做的好处是：

- 调试数据先在内存结构里稳定下来
- 后续就算不用 probe app，也能在别的 app 里复用这些统计

### 3.3 尽量不改 Chase core / decode_tile_with_plan 的语义

这次实验不建议动的地方：

- Chase 候选生成
- BCH 硬解码算法
- MUX 调度规则本身
- A/B merged 顺序本身

建议只做：

- 统计透出
- 结构扩展
- CSV 导出

这样更容易保证：

- 调试功能本身不会反过来影响 BER 行为

---

## 4. 需要改哪些代码

下面按“必须修改”和“建议新增”分组说明。

## 4.1 必须修改：`include/newcode/two_stream_shared_runner.hpp`

文件：

- [two_stream_shared_runner.hpp](/home/zsr71/projects/newcode_two_stream_shared_chase/include/newcode/two_stream_shared_runner.hpp:1)

### 目标

扩展 shared runner 的对外结果结构，让它不仅返回：

- quantization 统计
- produced/failed 汇总

还能够返回：

- tile 级样本
- hybrid 分类汇总
- row 级 map

### 计划修改内容

新增几组结构体。

第一组：tile 级样本结构

- `SharedTileSample`

建议字段包括：

- `invocation`
- `tile_index`
- `stream_rows_a`
- `stream_rows_b`
- `rows_total`
- `rows_early_stop`
- `rows_not_early_stop`
- `rows_hard_finish`
- `rows_need_siso_before_mux`
- `rows_soft_scheduled`
- `rows_unscheduled`
- `produced_rows`
- `failed_rows`
- `produced_rows_a`
- `produced_rows_b`
- `failed_rows_a`
- `failed_rows_b`

第二组：hybrid 分类统计结构

- `SharedHybridClassCount`

建议字段包括：

- `invocation`
- `tile_index`
- `rows_seen_by_hybrid`
- `class_none_count`
- `class_bch_hard_decoded_count`
- `class_clean_count`
- `class_parity_only_count`
- `class_one_main_count`
- `class_one_main_plus_parity_count`
- `class_two_main_count`
- `class_suspicious_count`
- `class_hard_fail_count`
- `deferred_candidate_count`
- `deferred_priority_0_count`
- `deferred_priority_1_count`
- `deferred_priority_2_count`
- `deferred_reclaimed_to_hard_finish_count`

第三组：row 级 map 结构

- `SharedRowMapEntry`

建议字段包括：

- `invocation`
- `tile_index`
- `merged_row`
- `stream_id`
- `source_local_row`
- `source_global_row`
- `early_stop_hit`
- `hybrid_class`
- `final_tag`
- `scheduled_for_soft`
- `produced_row`

然后把这些结构挂到：

- `Observability`

里，例如：

- `std::vector<SharedTileSample> tile_samples;`
- `std::vector<SharedHybridClassCount> hybrid_class_counts;`
- `std::vector<SharedRowMapEntry> row_map;`

### 为什么必须改这个文件

因为这一步决定了：

- shared runner 的调试结果是否成为正式 API

如果不先把结构定义在头文件里，后面 probe app 就只能重新钻进内部实现里临时拿数据，不利于长期维护。

---

## 4.2 必须修改：`src/two_stream_shared/two_stream_shared_runner.cpp`

文件：

- [two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:1)

### 目标

把当前 shared tile 处理过程中已经拿到的信息，真正填充进新定义的 observability 结构。

### 计划修改内容

#### 4.2.1 扩展 `SharedTileResult`

当前 `SharedTileResult` 已经有数量摘要，但还不够支撑三份 CSV。

需要为它补两类信息：

1. 当前 tile 的 row 级 map 数据
2. 当前 tile 的 hybrid class 汇总结果

如果不想让 `SharedTileResult` 太重，也可以：

- 新增一个更上层的 debug 结构
- 在 `run_shared_tile(...)` 内部组装后直接返回

但从代码改动量看，直接扩展 `SharedTileResult` 更顺手。

#### 4.2.2 在 `run_shared_tile(...)` 中收集 tile 样本

在
[two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:474)
之后，当前已经有：

- `dispatch_plan`
- `rows_need_siso_before_mux`
- `decoder_res`
- `shared_prep.slices`

这里正是组装三类调试信息的最佳位置。

要补的事情包括：

- 计算 `rows_not_early_stop`
- 读取 `dispatch_plan.rows_soft_scheduled`
- 读取 `dispatch_plan.rows_soft_unscheduled`
- 生成一条 `SharedTileSample`

#### 4.2.3 在 `run_shared_tile(...)` 中生成 row map

这里需要遍历 merged row 域，把下面这些信息逐行拼出来：

- `merged_row`
- 当前 row 属于 A 还是 B
- 在各自 slice 中的 local row 是多少
- 在原始 frame 里的 global row 是多少
- `early_stop_hit`
- `hybrid_class`
- `final_tag`
- `scheduled_for_soft`
- `produced_row`

这些信息的来源分别是：

- 来源流 / local row：来自 `shared_prep.slices`
- global row：来自 `shared_prep.merged.row_global_lookup`
- early-stop / hybrid / final tag / scheduled：来自 `dispatch_plan.rows[...]`
- produced：来自 `decoder_res.produced_rows[...]`

#### 4.2.4 在 `run_shared_tile(...)` 中生成 hybrid class count

这一步可以直接扫描：

- `dispatch_plan.hybrid_classes`

按枚举值做计数。

额外需要计的 deferred 相关列，当前还没有正式保存，见下一节。

#### 4.2.5 在 `decode_two_stream_shared_llr(...)` 中把每个 tile 的调试结果聚合到总 observability

在
[two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp:663)
附近，当前只把：

- `shared_core_stats`

往上累计。

这里需要再补：

- `tile_samples.push_back(...)`
- `hybrid_class_counts.push_back(...)`
- `row_map.insert(...)`

同时维护：

- `invocation` 计数器

建议这个计数器与当前 window/tile 主循环同生命周期，按每次 `run_shared_tile(...)` 调用递增。

### 为什么必须改这个文件

因为这是真正的 shared 主流程实现所在位置。

如果不在这里收集调试信息，后续 app 层很难在不重复执行逻辑的情况下准确恢复：

- row 去向
- hybrid 分类
- shared 预算竞争结果

---

## 4.3 必须修改：`src/rx/ofec/detail/ofec_tile_dispatch.ipp`

文件：

- [ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:1)

### 目标

把当前 hybrid prepass 里“只在局部变量里存在”的 deferred/backfill 统计正式挂进 `TileDispatchPlan`。

### 当前缺口

当前代码里：

- `hybrid_classes` 已经保留了
- `rows_hard_finish` 已经保留了
- `rows_soft_candidate / rows_soft_scheduled / rows_soft_unscheduled` 也已经保留了

但是这些信息还不够支撑第二份 CSV 里的 deferred 相关列。

因为在
[ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/rx/ofec/detail/ofec_tile_dispatch.ipp:175)
里，下面这些量目前只是局部变量：

- `deferred_candidates_total`
- `deferred_candidates_by_priority`
- `reclaimed`

等函数返回之后，这些信息就丢了。

### 计划修改内容

给 `TileDispatchPlan` 增加几项统计字段：

- `deferred_candidate_count`
- `deferred_priority_0_count`
- `deferred_priority_1_count`
- `deferred_priority_2_count`
- `deferred_reclaimed_to_hard_finish_count`

然后在 `run_hybrid_prepass(...)` 里填这些值。

### 为什么必须改这个文件

如果不改这里，就无法稳定回答：

- `TwoMain` 里有多少被保留成 deferred
- 最后又有多少被 reclaim 回 `HardFinish`

而这是第二份 CSV 非常关键的一组解释信息。

---

## 4.4 建议新增：`apps/ofec_two_stream_ber_window_probe.cpp`

文件：

- 新增 [ofec_two_stream_ber_window_probe.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_ber_window_probe.cpp)

### 目标

做一个专门的 two-stream probe 入口，风格参考：

- [ofec_ber_window_probe.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_ber_window_probe.cpp:1)

但执行对象改成：

- `newcode::two_stream_shared::run_two_stream_shared(...)`

### 为什么推荐新增而不是直接塞进现有 app

原因有三个：

1. 单流 probe 的结果结构是按单流 `PipelineResult` 设计的
2. 双流 shared probe 需要同时落盘 shared 内部三份表
3. one-shot app 更适合打印摘要，不适合承担多 seed / 多 `Eb/N0` / 多 CSV 管理

### 计划实现内容

这个 app 的职责应包括：

- 组织 `Eb/N0` 列表
- 组织 seed 列表
- 给 A/B 生成独立前端 seed
- 调用 `run_two_stream_shared(...)`
- 从 `result.observability` 里取出三类调试结果
- 写出三份 CSV
- 视需要再补：
  - A/B 的 per-window BER
  - A/B 的 per-tile BER
  - 日志摘要

### 参数设计建议

参数设计应参考：

- [ofec_two_stream_ber_window_probe参数说明.md](/home/zsr71/projects/newcode_two_stream_shared_chase/doc/ofec_two_stream_ber_window_probe参数说明.md:1)

第一版建议固定：

- 先只跑很少的 `Eb/N0`
- 很少的 seed
- 重点验证 CSV 结构和 shared 内部统计是否正确

---

## 4.5 建议修改：`CMakeLists.txt`

文件：

- [CMakeLists.txt](/home/zsr71/projects/newcode_two_stream_shared_chase/CMakeLists.txt:207)

### 目标

把新的 probe app 加进编译系统。

### 计划修改内容

仿照现有：

- `ofec_two_stream_shared_decoder`
- `ofec_two_stream_shared_sweep`

新增：

- `ofec_two_stream_ber_window_probe`

建议链接：

- `newcode_frontend`
- `decoder_plain`
- `decoder_ebchPF`

并复用：

- `src/ofec_single/ofec_single_dual_writer.cpp`
- `src/ofec_single/ofec_single_params.cpp`
- `src/two_stream_shared/shared_bch_hard_decode64.cpp`
- `src/two_stream_shared/two_stream_shared_runner.cpp`

如果 probe 里需要复用单流 probe 的部分辅助逻辑，也可以再考虑抽公共 helper，但第一版不建议为此大动结构。

---

## 4.6 可选修改：`apps/ofec_two_stream_shared_decoder.cpp`

文件：

- [ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp:1)

### 目标

不是让它写 CSV，而是：

- 临时加一点摘要打印
- 方便在 probe app 没写完前，快速 sanity check observability 结构

### 建议修改内容

在拿到 `result.observability` 后，可选打印：

- tile sample 总数
- hybrid class count 总数
- row map 总行数

以及若干聚合摘要，例如：

- total rows early-stop
- total rows hard-finish
- total rows unscheduled

### 为什么是可选

因为这部分不是最终实验交付物。

它主要是开发阶段用于：

- 快速看结构是不是有数据

---

## 5. 不建议改哪些代码

为了控制风险，这一轮不建议动下面这些地方：

### 5.1 不建议改 Chase 核心实现

例如：

- `chase256_plain_impl.ipp`
- `chase256_topk_pruned_impl.ipp`
- `chase256_global_pair_impl.ipp`

原因：

- 这次目标不是算法修正，而是调试可视化

### 5.2 不建议改 early-stop 判定逻辑

例如：

- `tile_early_stop_stats.ipp`
- `tile_early_stop_stats2.ipp`

原因：

- 这些逻辑当前已经能提供 `row_passed_flags`
- 足够做调试统计

### 5.3 不建议一开始就加 Chase candidate 级 CSV

原因：

- 当前三份 shared CSV 的粒度已经覆盖：
  - tile 总体
  - hybrid 分类
  - row 去向
- 如果再把 Chase candidate 明细混进来，第一版复杂度会明显上升

---

## 6. 推荐的分阶段实施顺序

## 6.1 第一阶段：先补 shared runner 结果结构

涉及文件：

- `include/newcode/two_stream_shared_runner.hpp`
- `src/two_stream_shared/two_stream_shared_runner.cpp`
- `src/rx/ofec/detail/ofec_tile_dispatch.ipp`

这一阶段完成后，目标是：

- `run_two_stream_shared(...)` 已经能在内存里返回三类调试结果
- 即使还没有 probe app，也可以在 one-shot app 中打印统计摘要

### 这一阶段验证方法

做一个固定参数的 one-shot run，检查：

- `tile_samples.size() > 0`
- `hybrid_class_counts.size() == tile_samples.size()`
- `row_map.size()` 等于所有 invocation 的 `rows_total` 之和

同时核对数量关系：

- `rows_total = stream_rows_a + stream_rows_b`
- `rows_not_early_stop = rows_total - rows_early_stop`
- `rows_soft_scheduled + rows_unscheduled = rows_need_siso_before_mux`
- `produced_rows + failed_rows = rows_total`

## 6.2 第二阶段：新增 two-stream probe app

涉及文件：

- `apps/ofec_two_stream_ber_window_probe.cpp`
- `CMakeLists.txt`

这一阶段完成后，目标是：

- 可以直接跑出三份 CSV
- 每个 invocation 都有稳定记录

### 这一阶段验证方法

用非常小的实验配置：

- 1 个 `Eb/N0`
- 1 个 seed
- 1 个 chunk

检查：

- 三份 CSV 都能生成
- 三份 CSV 的 `invocation / tile_index` 能对齐
- 第一份和第二份是一对一
- 第三份每个 invocation 的行数等于第一份里的 `rows_total`

## 6.3 第三阶段：扩大 seed / `Eb/N0` 范围，做真实诊断实验

这一阶段不再新增结构，而是用已经稳定的 probe：

- 跑多个 seed
- 跑 1 到 3 个 `Eb/N0`
- 对 A/B 交换实验做时域级对比

### 这一阶段重点看什么

- A/B 在 `row_map` 层面是否存在系统性偏差
- `TwoMain` / `HardFail` 的占比是否随 `Eb/N0` 改变
- 某些 tile 是否长期 `rows_unscheduled` 偏高

---

## 7. 每个文件的具体改动摘要

为了方便实际落地，下面给一个更简洁的“文件 -> 改动点”索引。

### 7.1 `include/newcode/two_stream_shared_runner.hpp`

新增：

- `SharedTileSample`
- `SharedHybridClassCount`
- `SharedRowMapEntry`

扩展：

- `Observability`

### 7.2 `src/two_stream_shared/two_stream_shared_runner.cpp`

扩展：

- `SharedTileResult`

新增逻辑：

- 在 `run_shared_tile(...)` 里生成 tile sample
- 在 `run_shared_tile(...)` 里生成 hybrid count
- 在 `run_shared_tile(...)` 里生成 row map
- 在 `decode_two_stream_shared_llr(...)` 里汇总这些结果

### 7.3 `src/rx/ofec/detail/ofec_tile_dispatch.ipp`

扩展：

- `TileDispatchPlan`

新增统计：

- deferred candidate 相关计数
- reclaimed 相关计数

### 7.4 `apps/ofec_two_stream_ber_window_probe.cpp`

新增：

- two-stream probe 主入口
- 三份 CSV 的表头和写文件逻辑
- seed / `Eb/N0` 调度

### 7.5 `CMakeLists.txt`

新增：

- `ofec_two_stream_ber_window_probe` 可执行目标

### 7.6 `apps/ofec_two_stream_shared_decoder.cpp`

可选扩展：

- 打印 observability 规模和若干摘要数

---

## 8. 风险点与对应处理

## 8.1 风险：row map 太大

第三份 CSV 可能会很大，尤其在：

- 多 seed
- 多 chunk
- 多 `Eb/N0`

时增长很快。

### 处理建议

第一版先限制：

- `Eb/N0` 点少
- seed 少
- chunk 少

必要时给 probe 加开关：

- 是否导出 row map

## 8.2 风险：调试结构改变运行时内存占用

如果直接把所有 row map 都长期存到 `Observability`，内存会变大。

### 处理建议

第一版可以接受这点代价，因为 probe 本来就是调试工具。

如果后续规模扩大，再考虑：

- probe 模式才开启 row map
- one-shot / sweep 默认不保存 row map

## 8.3 风险：数量关系不一致

例如：

- `rows_soft_scheduled + rows_unscheduled != rows_need_siso_before_mux`

这通常意味着：

- 统计口径拿错位置
- 或某些 tag 变更没有同步更新计数

### 处理建议

第一阶段就加断言式自检或 debug 校验，优先在内存里把关系对齐，再写 CSV。

---

## 9. 最小可交付版本

如果按“最小能开始分析 shared 内部行为”的标准，建议第一批必须交付的是：

1. `SharedTileSample`
2. `SharedRowMapEntry`
3. `apps/ofec_two_stream_ber_window_probe.cpp`

第二份 `SharedHybridClassCount` 虽然也重要，但如果要进一步压缩第一批复杂度，也可以稍晚一点补。

不过从当前代码基础看，既然：

- `hybrid_classes` 已经有了
- deferred 计数只差几处结构补充

所以实际更推荐一次把三份表都做齐。

---

## 10. 一句话结论

这次实验如果按最小侵入路线推进，真正需要动的核心代码只有三块：

1. `two_stream_shared_runner.hpp`
2. `two_stream_shared_runner.cpp`
3. `ofec_tile_dispatch.ipp`

然后再新增一块独立实验入口：

4. `apps/ofec_two_stream_ber_window_probe.cpp`

以及一处编译接线：

5. `CMakeLists.txt`

整个思路不是“改 shared 算法”，而是：

> 先把 shared 内部调度与分类信息变成正式结构，再用独立 probe app 把它们稳定导出成三份 CSV。  

这样后续无论是定位 A/B 不对称，还是分析 hybrid / MUX 行为，都会有一个足够干净、足够稳定的观测面。
