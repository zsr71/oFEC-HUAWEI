# `two-stream shared` 调试 CSV 设计说明

本文说明当前已经落地的一组 `two-stream shared` 调试 CSV。

这些 CSV 的目标不是替代最终的 BER 结果表，而是专门服务于：

- 分析 shared 64-row 域内部到底发生了什么
- 解释 A/B 两路为什么会出现不对称
- 观察 early-stop、hybrid、MUX、shared soft decode 在每次 tile 调用里的分流结果
- 让后续 Matlab / Python 分析可以直接基于结构化表格做筛选、聚合和可视化

本文描述的三份 CSV 分别是：

1. `per_invocation_shared_tile_samples.csv`
2. `per_invocation_shared_hybrid_class_counts.csv`
3. `per_invocation_shared_row_map.csv`

这三份表的关系是：

- 第一份看“这一整个 shared tile 调用的总体摘要”
- 第二份看“这个 shared tile 里 hybrid 分类具体分成了多少种情况”
- 第三份看“这个 shared tile 里每一行最终到底被分到了哪里”

可以把它们理解成三层视角：

- tile 级摘要
- tile 级分类细分
- row 级明细

---

## 1. 总体设计原则

### 1.1 这三份 CSV 都面向“单次 shared tile 调用”

这里的 `invocation` 指的是：

- 在一次具体运行中
- 某个 seed / chunk / `Eb/N0` 条件下
- shared runner 每进入一次 merged 64-row tile 路径
- 就记一条 invocation

也就是说，本文的“每次进来”对应的就是：

- 一次 `merge_prepared_tiles(...)`
- 一次 `detect_tile_early_stop(...)`
- 一次 `build_tile_dispatch_plan(...)`
- 一次 `run_mux_on_soft_candidates(...)`
- 一次 `decode_tile_with_plan(...)`

完成后的整套 shared tile 处理

### 1.2 三份表之间要能直接对齐

因此三份 CSV 都应共享一组最基本的键列：

- `run_id`
- `ebn0_db`
- `seed_index`
- `chunk_index`
- `invocation`
- `tile_index`

含义是：

- `run_id` 用来区分不同批次的实验运行
- `ebn0_db` 用来区分不同噪声点
- `seed_index` 用来区分同一个 `Eb/N0` 下的不同 seed
- `chunk_index` 用来区分 low-BER sweep 里的不同 chunk
- `invocation` 用来区分一次运行过程中第几次 shared tile 调用
- `tile_index` 用来区分当前 window 内的第几个 tile

如果是 one-shot app 而不是 sweep，也仍然建议保留：

- `chunk_index`

只是固定写成 `0`。

### 1.3 先做“足够稳定的结构化信息”，不急着一次做完所有 debug 细节

第一版重点应放在：

- 让 shared 内部关键数量关系可见
- 让 A/B 两路的 row 去向可回放
- 让 hybrid 分类细分可统计

而不是一开始就把 Chase 每个 candidate 的所有细节都并到同一个大 CSV 里。

Chase candidate 明细如果后续要看，更适合单独做第四份表，而不是塞进这里三份里。

### 1.4 四段流程与采集时序

为了正确理解每一列，先把一次 shared tile invocation 的内部时序固定下来。

当前实现可以近似分成下面四段：

1. `early-stop` 阶段
2. `hybrid` 分类 / hard-finish / deferred 回收阶段
3. `soft candidate` 保留阶段
4. `MUX` 裁剪 + 最终 decode 阶段

更贴近代码的顺序是：

1. `merge_prepared_tiles(...)`
2. `detect_tile_early_stop(...)`
3. `build_tile_dispatch_plan(...)`
4. `run_hybrid_prepass(...)`
5. `count_soft_candidates_before_mux(...)`
6. `run_mux_on_soft_candidates(...)`
7. `decode_tile_with_plan(...)`

这四段和代码的对应关系如下：

- 第 1 段 `early-stop`
  - 对应 `detect_tile_early_stop(...)`
  - 输出“哪些 row 已经 early-stop”
  - 这一步结束后，这些 row 已经被固定成 `EarlyStopAction`

- 第 2 段 `hybrid 分类 / hard-finish / deferred 回收`
  - 对应 `build_tile_dispatch_plan(...)` 里的 `run_hybrid_prepass(...)`
  - 只处理“没有 early-stop、且当前 tag 仍然是 `SoftDecode` 的 row”
  - 这一步内部既做分类，也做 deferred priority 统计，还会在条件满足时把一部分 deferred row 收回成 `HardFinish`

- 第 3 段 `soft candidate 保留`
  - 对应 `count_soft_candidates_before_mux(...)`
  - 统计的是“经过 early-stop 和 hybrid 之后，还剩多少 row 仍然保留在 soft 候选池里”
  - 这一步是 MUX 之前的快照，不包含 MUX 裁剪结果

- 第 4 段 `MUX 裁剪 + 最终 decode`
  - 对应 `run_mux_on_soft_candidates(...)` 和 `decode_tile_with_plan(...)`
  - `run_mux_on_soft_candidates(...)` 决定哪些 soft row 真正保留为 `SoftDecode`，哪些变成 `Unscheduled`
  - `decode_tile_with_plan(...)` 再基于最终 tag 真正产出 `produced_row`

因此，下面三份 CSV 里看到的后半段列，必须区分它们是：

- 第 1 段结束后的结果
- 第 2 段结束后的结果
- 第 3 段，也就是 `MUX` 之前的快照
- 第 4 段，也就是 `MUX` 之后 / decode 之后的最终结果

---

## 2. 三份 CSV 的整体关系

推荐把三份表理解成：

### 2.1 `per_invocation_shared_tile_samples.csv`

这是 tile 级总表。

一行代表：

- 一个具体 invocation
- 在某个 shared merged tile 上
- 整体发生了什么

这个表最适合回答的问题是：

- 这一轮 shared tile 一共多少行
- 有多少行 early-stop
- 有多少行变成 hard-finish
- 有多少行进入 soft 候选
- 最后有多少行被 MUX 裁掉
- produced / failed 的结果如何
- A/B 两路分别产出了多少 / 失败了多少

### 2.2 `per_invocation_shared_hybrid_class_counts.csv`

这是 hybrid 分类细分表。

一行同样代表：

- 一个具体 invocation
- 在某个 shared merged tile 上
- hybrid prepass 的分类统计结果

这个表最适合回答的问题是：

- hybrid 里到底是 `Clean` 多，还是 `OneMain` 多
- `TwoMain` 的数量有多少
- `HardFail` 的数量有多少
- 哪些分类会进入 deferred/backfill 候选
- 最后真正被 reclaim 成 hard-finish 的数量有多少

### 2.3 `per_invocation_shared_row_map.csv`

这是 row 级去向明细表。

一行代表：

- 一个具体 merged row
- 它来自 A 还是 B
- 它在 shared path 里是否 early-stop
- hybrid 给它打了什么分类
- 它最终走的是 `EarlyStopAction`、`HardFinish`、`SoftDecode` 还是 `Unscheduled`
- 它有没有抢到 soft decode 预算

这个表最适合回答的问题是：

- 为什么 A/B 交换后结果不对称
- 某一类 row 是否更容易被裁掉
- 某些 tile 是否对 A 或 B 更偏
- `TwoMain` 候选最后是走 hard-finish 还是被保留去抢 SISO

---

## 3. 第一份 CSV：`per_invocation_shared_tile_samples.csv`

## 3.1 一行代表什么

每一行代表：

- 一次 shared tile 调用的整体摘要

换句话说，这一行不是单个 row 的信息，而是整个 merged tile 的总计信息。

如果当前 shared tile 是：

- 32 行来自 A
- 32 行来自 B

那么这一行统计的是这 64 行整体的分流与结果。

## 3.2 推荐列定义

推荐列如下：

```text
run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,
stream_rows_a,stream_rows_b,rows_total,
rows_early_stop,rows_not_early_stop,
rows_hard_finish,rows_need_siso_before_mux,
rows_soft_scheduled,rows_unscheduled,
produced_rows,failed_rows,
produced_rows_a,produced_rows_b,
failed_rows_a,failed_rows_b
```

## 3.3 各列含义

### 标识列

`run_id`

- 当前整次实验运行的编号
- 同一次运行中的三份 CSV 都应相同

`ebn0_db`

- 当前 invocation 所属的 `Eb/N0`

`seed_index`

- 当前 invocation 所属的 seed 编号

`chunk_index`

- 当前 invocation 所属的 chunk 编号
- one-shot 模式下可固定为 `0`

`invocation`

- 当前运行中的第几次 shared tile 调用
- 建议从 `0` 或 `1` 单调递增

`tile_index`

- 当前 window 内的 tile 编号

### shared 规模列

`stream_rows_a`

- 当前 merged tile 中来自 Stream A 的 row 数

`stream_rows_b`

- 当前 merged tile 中来自 Stream B 的 row 数

`rows_total`

- 当前 merged tile 的总 row 数
- 通常应满足：
  `rows_total = stream_rows_a + stream_rows_b`

### early-stop / hybrid / MUX 摘要列

`rows_early_stop`

- 命中 early-stop 的 row 数
- 这些行会直接走 `EarlyStopAction`
- 采集时刻：
  - 第 1 段结束后就已经确定
  - 后续 hybrid / MUX / decode 都不会再改变这个数量

`rows_not_early_stop`

- 未命中 early-stop 的 row 数
- 通常应满足：
  `rows_not_early_stop = rows_total - rows_early_stop`
- 采集时刻：
  - 第 1 段结束后就已经确定
  - 它也是第 2 段 hybrid 的输入规模

`rows_hard_finish`

- 在 hybrid prepass 后直接变成 `HardFinish` 的 row 数
- 这里包含两部分：
  - 分类后立刻确定为 `HardFinish` 的 row
  - 先进入 deferred，随后又被回收成 `HardFinish` 的 row
- 采集时刻：
  - 第 2 段结束后统计
  - 也就是 hybrid 分类和 deferred reclaim 都做完之后、MUX 开始之前

`rows_need_siso_before_mux`

- 经过 early-stop 和 hybrid prepass 之后，仍然需要去竞争 SISO 的 row 数
- 也就是 MUX 之前的 soft candidate 数
- 采集时刻：
  - 第 3 段统计
  - 更精确地说，是第 2 段刚结束、`run_mux_on_soft_candidates(...)` 还没开始时的快照
- 这个值不是最终会进入 soft decode 的数量，而是“进入 MUX 竞争池之前”的数量

`rows_soft_scheduled`

- 最终抢到 soft decode 预算、会进入 soft path 的 row 数
- 采集时刻：
  - 第 4 段中 `run_mux_on_soft_candidates(...)` 执行完成后统计
  - 它表示 MUX 裁剪之后仍保留 `SoftDecode` tag 的 row 数

`rows_unscheduled`

- 经过 MUX 后被裁掉、不会进入本轮 soft decode 的 row 数
- 采集时刻：
  - 第 4 段中 `run_mux_on_soft_candidates(...)` 执行完成后统计
  - 这些 row 在这一步会被改写成 `Unscheduled`

通常这些数量关系应大致满足：

```text
rows_soft_scheduled + rows_unscheduled = rows_need_siso_before_mux
```

### shared decode 结果列

`produced_rows`

- 当前 invocation 最终成功产出 decoder 输出的 row 数
- 采集时刻：
  - 第 4 段末尾，也就是 `decode_tile_with_plan(...)` 真正跑完之后统计

`failed_rows`

- 当前 invocation 最终没有产出 decoder 输出的 row 数
- 采集时刻：
  - 第 4 段末尾，与 `produced_rows` 同时统计

`produced_rows_a`

- 这些成功产出的 row 中，来自 Stream A 的数量
- 采集时刻：
  - 第 4 段末尾，在 decode 结果按 A/B slice 写回时统计

`produced_rows_b`

- 这些成功产出的 row 中，来自 Stream B 的数量
- 采集时刻：
  - 第 4 段末尾，在 decode 结果按 A/B slice 写回时统计

`failed_rows_a`

- 失败 row 中，来自 Stream A 的数量
- 采集时刻：
  - 第 4 段末尾，在 decode 结果按 A/B slice 写回时统计

`failed_rows_b`

- 失败 row 中，来自 Stream B 的数量
- 采集时刻：
  - 第 4 段末尾，在 decode 结果按 A/B slice 写回时统计

通常应满足：

```text
produced_rows + failed_rows = rows_total
produced_rows_a + produced_rows_b = produced_rows
failed_rows_a + failed_rows_b = failed_rows
```

## 3.4 适合看什么

这张表适合直接看：

- 某个 tile 是否 early-stop 特别多
- 某个 tile 是否 soft 候选很多但大部分被 unscheduled
- A/B 两路 produced/failed 是否平衡
- 共享预算收紧后，`rows_unscheduled` 是否明显上升

---

## 4. 第二份 CSV：`per_invocation_shared_hybrid_class_counts.csv`

## 4.1 一行代表什么

每一行代表：

- 一个 shared tile invocation 内
- hybrid prepass 对所有“未 early-stop 且进入 prepass 的 row”
- 打出来的分类计数

它不关心单个 row 是哪一条，而关心：

- 这一整块 tile 里各种 hybrid class 各有多少

## 4.2 推荐列定义

推荐列如下：

```text
run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,
rows_seen_by_hybrid,
class_none_count,
class_bch_hard_decoded_count,
class_clean_count,
class_parity_only_count,
class_one_main_count,
class_one_main_plus_parity_count,
class_two_main_count,
class_suspicious_count,
class_hard_fail_count,
deferred_candidate_count,
deferred_priority_0_count,
deferred_priority_1_count,
deferred_priority_2_count,
deferred_reclaimed_to_hard_finish_count
```

## 4.3 各列含义

### 标识列

`run_id, ebn0_db, seed_index, chunk_index, invocation, tile_index`

- 与第一份表完全一致
- 用来直接 join

### hybrid 输入规模列

`rows_seen_by_hybrid`

- 真正进入 hybrid prepass 的 row 数
- 它不是 `rows_total`，因为：
  - early-stop 已经命中的 row 不会再进入 hybrid
- 采集时刻：
  - 第 2 段内部统计
  - 每当一行满足“未 early-stop 且当前仍是 `SoftDecode`”并真正进入 hybrid prepass 时，就会累计
- 它表示 hybrid 实际看到了多少行，而不是 tile 总行数

因此通常应满足：

```text
rows_seen_by_hybrid = rows_not_early_stop
```

### hybrid 分类列

`class_none_count`

- 理论上通常应为 `0`
- 仅当某些 row 没有被正确打 class 时才可能非 0
- 采集时刻：
  - 第 2 段结束后，根据每个进入 hybrid 的 row 最终写回的 `hybrid_class` 汇总

`class_bch_hard_decoded_count`

- LegacyHardDecode 模式下，直接调用传统 BCH 硬译码成功的数量
- 采集时刻：
  - 第 2 段结束后汇总
  - 属于 hybrid 分类的最终结果之一

`class_clean_count`

- 被判成已经是合法码字的 row 数
- 采集时刻：
  - 第 2 段结束后汇总

`class_parity_only_count`

- 只需要修 overall parity 的 row 数
- 采集时刻：
  - 第 2 段结束后汇总

`class_one_main_count`

- 被判成主体 1 错的 row 数
- 采集时刻：
  - 第 2 段结束后汇总

`class_one_main_plus_parity_count`

- 被判成主体 1 错 + parity 也要一起修的 row 数
- 采集时刻：
  - 第 2 段结束后汇总

`class_two_main_count`

- 被判成主体 2 错的 row 数
- 采集时刻：
  - 第 2 段结束后汇总

`class_suspicious_count`

- 预留给“分类可疑但不能直接 hard-finish”的 row 数
- 采集时刻：
  - 第 2 段结束后汇总

`class_hard_fail_count`

- 明确不满足当前 hard path 条件、只能继续走 soft path 的 row 数
- 采集时刻：
  - 第 2 段结束后汇总

### deferred / backfill 相关列

`deferred_candidate_count`

- 本来已经可以 hard-finish，但由于 backfill 策略，先被保留为 soft 候选的 row 数
- 采集时刻：
  - 第 2 段内部统计
  - 某行在 hybrid 分类时已经满足 hard-finish，但被判定为可以先 defer 给 soft path 时立即累计
- 这发生在 MUX 之前，属于 hybrid 内部的“先保留”动作

`deferred_priority_0_count`

- deferred 候选里优先级 0 的数量
- 采集时刻：
  - 第 2 段内部统计
  - 在某个 deferred row 被放入 priority 0 队列时累计

`deferred_priority_1_count`

- deferred 候选里优先级 1 的数量
- 采集时刻：
  - 第 2 段内部统计
  - 在某个 deferred row 被放入 priority 1 队列时累计

`deferred_priority_2_count`

- deferred 候选里优先级 2 的数量
- 采集时刻：
  - 第 2 段内部统计
  - 在某个 deferred row 被放入 priority 2 队列时累计

`deferred_reclaimed_to_hard_finish_count`

- 最终从 deferred 候选中又被收回、真正变成 `HardFinish` 的数量
- 采集时刻：
  - 第 2 段末尾统计
  - 等所有 deferred 候选都收集完、再结合 SISO 预算做一次 reclaim 后得到
- 它发生在 `rows_need_siso_before_mux` 统计之前

## 4.4 适合看什么

这张表适合直接分析：

- `TwoMain` 在不同 `Eb/N0` 下的占比
- `HardFail` 是否过高
- 开启 `OneAndTwoErrorPriority` 后，deferred 候选是否显著增多
- 某些 tile 的 `Clean / ParityOnly / OneMain / TwoMain` 结构是否和别的 tile 很不一样

---

## 5. 第三份 CSV：`per_invocation_shared_row_map.csv`

## 5.1 一行代表什么

每一行代表：

- 一个具体 merged row

它是三份表里最细的一份。

如果某个 invocation 有 64 个 merged row，那么这张表就会有 64 行对应它。

## 5.2 推荐列定义

推荐列如下：

```text
run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,
merged_row,stream_id,stream_label,source_local_row,source_global_row,
early_stop_hit,hybrid_class,final_tag,scheduled_for_soft,
produced_row
```

如果后续需要更细，也可以继续扩展：

```text
failed_row,
deferred_candidate,deferred_priority,
soft_candidate_before_mux,
soft_scheduled_after_mux
```

## 5.3 各列含义

### 标识列

`run_id, ebn0_db, seed_index, chunk_index, invocation, tile_index`

- 与前两张表完全一致

### row 来源列

`merged_row`

- 当前 row 在 merged tile 里的行号

`stream_id`

- 当前 row 来自哪一路
- 建议约定：
  - `0` = A
  - `1` = B

`stream_label`

- 当前 row 的来源流标签
- 建议写成：
  - `A`
  - `B`

`source_local_row`

- 当前 row 在各自 slice 内部的局部行号

`source_global_row`

- 当前 row 回到原始 frame / window 语义下的全局行号

### 调度语义列

`early_stop_hit`

- `1` 表示命中 early-stop
- `0` 表示未命中
- 采集时刻：
  - 第 1 段结束后写入 row map
  - 后续阶段不会再改这个布尔值

`hybrid_class`

- hybrid prepass 给这行打的分类标签
- 建议使用字符串枚举而不是数字，便于直接读：
  - `none`
  - `bch_hard_decoded`
  - `clean`
  - `parity_only`
  - `one_main`
  - `one_main_plus_parity`
  - `two_main`
  - `suspicious`
  - `hard_fail`
- 采集时刻：
  - 第 2 段结束后写入
  - 对于命中 early-stop 的 row，这里通常保持 `none`，因为它根本没有进入 hybrid
  - 对于未 early-stop 的 row，这里反映的是 hybrid 最终写回的分类结果

`final_tag`

- 这行最终走到哪条执行路径
- 建议值为：
  - `early_stop_action`
  - `hard_finish`
  - `soft_decode`
  - `unscheduled`
- 采集时刻：
  - 第 4 段开始前后都可能变化，但 row map 里记录的是最终值
  - 也就是 early-stop、hybrid、MUX 全部完成之后的最终调度结果
- 它的来源阶段是：
  - 第 1 段可能把 row 定成 `early_stop_action`
  - 第 2 段可能把 row 定成 `hard_finish`
  - 第 4 段的 MUX 可能把保留下来的 row 维持为 `soft_decode`，也可能改成 `unscheduled`

`scheduled_for_soft`

- `1` 表示最终抢到了 soft decode 预算
- `0` 表示没有
- 采集时刻：
  - 第 4 段中 `run_mux_on_soft_candidates(...)` 完成后写入
  - 它只反映 MUX 之后是否拿到了 soft 预算，不反映 MUX 之前是否曾经是 soft candidate

这里通常有以下关系：

- `final_tag = soft_decode` 时，`scheduled_for_soft` 应为 `1`
- `final_tag = unscheduled` 时，`scheduled_for_soft` 应为 `0`

### 最终结果列

`produced_row`

- `1` 表示这行最终有 decoder 输出
- `0` 表示没有
- 采集时刻：
  - 第 4 段末尾，在 `decode_tile_with_plan(...)` 真正跑完以后写入
  - 这是最靠后的“最终结果列”

如果后续你们更喜欢显式列，也可以增加：

`failed_row`

- 直接写成 `1 - produced_row`

## 5.4 适合看什么

这张表最适合做下面这些分析：

- A/B 两路 row 在同一 tile 里是否存在系统性偏差
- `TwoMain` 行在 A 路和 B 路的命运是否不同
- 某些 merged row 位置是否更容易被 MUX 裁掉
- 交换 A/B 后，是否只是 `stream_label` 交换了，但 `final_tag` 分布变了

---

## 6. 三份 CSV 之间怎么联合使用

推荐的分析顺序是：

### 第一步：先看第一份总表

先定位问题大概在哪些 invocation / tile：

- 哪些 invocation 的 `rows_unscheduled` 很高
- 哪些 invocation 的 `produced_rows_a` 和 `produced_rows_b` 差别很大

### 第二步：再看第二份 hybrid 分类表

确认这些异常 invocation 的 hybrid 结构是否特殊：

- 是不是 `TwoMain` 特别多
- 是不是 `HardFail` 特别多
- 是不是 deferred 候选太多

### 第三步：最后看第三份 row map

如果还要解释 A/B 不对称或某种分类的最终去向，再下钻到 row 级：

- 哪些 row 属于 A
- 哪些 row 属于 B
- 谁被 early-stop
- 谁被 hard-finish
- 谁去抢 soft budget
- 谁被 unscheduled

---

## 7. 一个完整例子

下面给一个简单例子。

假设某次运行里：

- `run_id = 20260604_143500`
- `ebn0_db = 3.05`
- `seed_index = 0`
- `chunk_index = 0`
- `invocation = 17`
- `tile_index = 3`

这个 invocation 里共有 64 行：

- A 路 32 行
- B 路 32 行

处理后出现以下情况：

- 18 行 early-stop
- 剩余 46 行进入 hybrid
- hybrid 里：
  - 8 行 `Clean`
  - 3 行 `ParityOnly`
  - 5 行 `OneMain`
  - 2 行 `OneMainPlusParity`
  - 4 行 `TwoMain`
  - 24 行 `HardFail`
- 其中 4 行 `TwoMain` 因 backfill 先保留为 deferred 候选
- 最后回收了其中 2 行做 `HardFinish`
- MUX 前一共还有 26 行需要抢 soft
- 预算只够 16 行
- 所以最终：
  - 16 行 `SoftDecode`
  - 10 行 `Unscheduled`
- 最后 produced 情况：
  - 总 produced = 61
  - 总 failed = 3
  - A produced = 30
  - B produced = 31

### 7.1 第一份表中的这一行

```csv
run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,stream_rows_a,stream_rows_b,rows_total,rows_early_stop,rows_not_early_stop,rows_hard_finish,rows_need_siso_before_mux,rows_soft_scheduled,rows_unscheduled,produced_rows,failed_rows,produced_rows_a,produced_rows_b,failed_rows_a,failed_rows_b
20260604_143500,3.05,0,0,17,3,32,32,64,18,46,20,26,16,10,61,3,30,31,2,1
```

这里的 `rows_hard_finish = 20` 表示：

- 不只是 `Clean/ParityOnly/OneMain/OneMainPlusParity`
- 还包括最终被 reclaim 成 `HardFinish` 的 deferred 候选

### 7.2 第二份表中的这一行

```csv
run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,rows_seen_by_hybrid,class_none_count,class_bch_hard_decoded_count,class_clean_count,class_parity_only_count,class_one_main_count,class_one_main_plus_parity_count,class_two_main_count,class_suspicious_count,class_hard_fail_count,deferred_candidate_count,deferred_priority_0_count,deferred_priority_1_count,deferred_priority_2_count,deferred_reclaimed_to_hard_finish_count
20260604_143500,3.05,0,0,17,3,46,0,0,8,3,5,2,4,0,24,4,0,0,4,2
```

这里表示：

- hybrid 一共看了 46 行
- 其中 24 行 `HardFail`
- 4 行 `TwoMain`
- 这 4 行全部先被保留成 deferred 候选
- 最后回收其中 2 行做 `HardFinish`

### 7.3 第三份表中的部分示例行

下面只列出 8 行示意，不展开全部 64 行：

```csv
run_id,ebn0_db,seed_index,chunk_index,invocation,tile_index,merged_row,stream_id,stream_label,source_local_row,source_global_row,early_stop_hit,hybrid_class,final_tag,scheduled_for_soft,produced_row
20260604_143500,3.05,0,0,17,3,0,0,A,0,384,1,none,early_stop_action,0,1
20260604_143500,3.05,0,0,17,3,1,1,B,0,384,0,clean,hard_finish,0,1
20260604_143500,3.05,0,0,17,3,2,0,A,1,385,0,parity_only,hard_finish,0,1
20260604_143500,3.05,0,0,17,3,3,1,B,1,385,0,one_main,hard_finish,0,1
20260604_143500,3.05,0,0,17,3,4,0,A,2,386,0,two_main,soft_decode,1,1
20260604_143500,3.05,0,0,17,3,5,1,B,2,386,0,two_main,unscheduled,0,0
20260604_143500,3.05,0,0,17,3,6,0,A,3,387,0,hard_fail,soft_decode,1,1
20260604_143500,3.05,0,0,17,3,7,1,B,3,387,0,hard_fail,unscheduled,0,0
```

从这几行就能直接看出：

- merged row 0 是 A 路，命中 early-stop
- merged row 1 是 B 路，被判成 `Clean`，直接 hard-finish
- merged row 4 是 A 路，`TwoMain`，最终抢到了 soft
- merged row 5 是 B 路，同样 `TwoMain`，但被 unscheduled

这正是第三份表最有价值的地方：

- 它能把“同一类 row 在 shared 资源竞争下最终命运不同”这件事直接摊开

---

## 8. 推荐的最小落地顺序

如果后面要实现，建议按下面顺序逐步落地：

1. 先做 `per_invocation_shared_tile_samples.csv`
2. 再做 `per_invocation_shared_row_map.csv`
3. 最后补 `per_invocation_shared_hybrid_class_counts.csv`

原因是：

- 第一份最容易直接利用现有计数器
- 第三份对定位 A/B 不对称最有帮助
- 第二份虽然也重要，但它依赖把 hybrid 分类计数正式汇总出来

---

## 9. 一句话总结

这三份 CSV 的职责分工应当是：

- `per_invocation_shared_tile_samples.csv`：看整块 tile 的数量关系
- `per_invocation_shared_hybrid_class_counts.csv`：看 hybrid 分类结构
- `per_invocation_shared_row_map.csv`：看每个 merged row 的最终去向

三份表一起用，才能真正把 `two-stream shared` 内部的：

- early-stop
- hybrid
- MUX
- shared soft decode
- A/B 不对称

这些现象从“感觉上有问题”变成“可以直接查表解释”。
