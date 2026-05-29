# 双虚拟流共享Chase方案验证分阶段实施建议

本文讨论的是：

> 如果希望基于当前仓库已有代码，验证“基于双虚拟流独立状态、共享 Chase 资源池”的架构方案是否有价值，应该分几个阶段推进，以及每个阶段如何实现。

这里的重点不是直接改 RTL，也不是一上来就在主链路里做两路独立实现，而是尽量复用当前软件仿真框架，先把关键问题分层验证清楚。

## 1. 为什么必须分阶段验证

两路独立 + 共享 Chase 方案的风险，不在于想法本身听起���不合理，而在于它同时叠加了多件事情：

- 双虚拟流划分
- 共享 Chase 资源池
- 全局调度
- 结果回写与重排序
- 拥塞退化策略

如果这些内容一次性全做进代码，最后即使看到：

- BER 变好了
- burst 变少了
- Chase 峰值变平了

也很难回答到底是哪一部分起了作用。

因此，建议把验证拆成 4 个阶段：

1. 先验证当前单流代码里，瞬时 Chase 负载是否真的存在 burst。
2. 再做离线“双流重放”，验证统计复用是否真的能削峰。
3. 再做“软件级共享资源池调度仿真”，验证队列、deadline miss 和拥塞行为。
4. 最后才考虑把双流结构以最小侵入方式接入现有主链路。

这样推进有两个好处：

- 每一阶段都能独立回答一个核心问题。
- 只有前一阶段成立，后一阶段才值得投入工程成本。

## 2. 当前代码里已经可以直接复用的观测点

当前仓库其实已经具备了不少和“共享 Chase 资源紧张”相关的统计基础，不需要从零开始。

### 2.1 tile 级资源相关统计

在 [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp:331) 到 [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp:447) 这段逻辑中，当前每个 tile 已经会统计：

- `rows_need_siso_before_mux`
- `rows_hard_finish`
- `rows_unscheduled`

这几个量的意义非常重要：

- `rows_need_siso_before_mux`
  - 当前 tile 在 MUX 之前，真正还想进入 soft path 的行数
  - 可以把它视作“当前 tile 的 Chase 请求需求”
- `rows_unscheduled`
  - 由于当前 SISO 预算不足，被裁掉没能进入 soft path 的行数
  - 可以把它视作“当前 tile 资源不足的直接症状”
- `rows_hard_finish`
  - 被 hybrid prepass 提前收走的行数
  - 它会影响后续真正进入 Chase 的需求规模

### 2.2 invocation 级样本输出

[ofec_single_runner.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_runner.cpp:7) 已经支持导出 invocation 级 tile 样本 CSV，目前字段是：

```text
invocation,tile_index,rows_total,rows_passed,
rows_hard_finish,rows_need_siso_before_mux,rows_unscheduled
```

这意味着：

- 我们已经能拿到按时间顺序排列的 tile 级资源需求样本。
- 这份 CSV 本质上就是验证双流统计复用最关键的输入之一。

### 2.3 hybrid / dispatch plan 已经提供了软路径候选定义

在 [ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_dispatch.ipp:1) 中，当前代码已经把 row 最终分成：

- `EarlyStopAction`
- `HardFinish`
- `SoftDecode`
- `Unscheduled`

这使得“哪些行真正构成 Chase 请求”在软件侧已经有了比较清楚的语义基础，不必重新定义。

## 3. 推荐的 4 个验证阶段

## 3.1 第一阶段：验证单流下是否真的存在瞬时 Chase burst

### 3.1.1 目标

先回答最基本的问题：

> 当前系统的 BER burst，是否真的和某些 tile / invocation 上的 Chase 请求峰值有关。

这一阶段不引入双流，不改主算法，只看单流统计。

### 3.1.2 建议指标

基于当前代码，先统计下面这些量的时间序列：

- `rows_need_siso_before_mux`
- `rows_unscheduled`
- `rows_hard_finish`
- `rows_total`
- `tile_index`
- `invocation`

然后重点看：

- `rows_need_siso_before_mux` 的峰值和长尾
- `rows_unscheduled` 是否以 burst 形式出现
- `rows_unscheduled > 0` 的 invocation 是否和 BER 恶化时段对齐
- 某些 tile 是否明显比其他 tile 更容易出现高负载

### 3.1.3 如何实现

当前代码基本已经能做：

1. 使用 `ofec_single` 或 `ofec_ber_window_probe` 跑固定配置。
2. 打开 tile 样本导出。
3. 收��� `data/early_stop_hist/...csv`。
4. 在 Python / Matlab 中按 `invocation` 画时序图和直方图。

如果想再补一点点代码，让验证更扎实，建议增加以下字段：

- `ebn0_db`
- `bitgen_seed`
- `channel_seed`
- `siso_active_for_tile`
- `decoder_name`

这样后面做多组 trace 汇总时不会丢上下文。

### 3.1.4 基于当前代码的推荐入口：优先使用 `ofec_ber_window_probe`

结合当前仓库现状，第一阶段更推荐优先使用
[apps/ofec_ber_window_probe.cpp](/home/zsr71/projects/newcode/apps/ofec_ber_window_probe.cpp:1)，
而不是直接从 `ofec_single` 开始。

原因不是名字里有 `probe`，而是它现在已经天然输出了更适合做“瞬时需求峰值分析”的几类 CSV：

- `per_seed_tile_early_stop_samples.csv`
- `per_seed_per_tile.csv`
- `per_seed_per_window.csv`
- `aggregated_tile_ber.csv`

其中最关键的是：

- `per_seed_tile_early_stop_samples.csv`

这份表当前已经按 `seed + invocation + tile_index` 粒度输出：

- `rows_total`
- `rows_passed`
- `rows_hard_finish`
- `rows_need_siso_before_mux`
- `rows_unscheduled`

对于第一阶段，建议把下面两个量分开看：

- `need_decode_after_early_stop = rows_total - rows_passed`
  - 它反映的是 early-stop 之后仍需要继续处理的行数
  - 更接近早期文档里“need_decode 峰值”的口径
- `rows_need_siso_before_mux`
  - 它反映的是 hybrid prepass 之后、MUX 之前，真正仍要去竞争 soft/SISO 资源的行数
  - 从当前代码语义看，它更接近“Chase 请求需求”的主指标

因此，第一阶段建议的主分析口径是：

1. 先用 `need_decode_after_early_stop` 看 early-stop 之后是否存在时域峰值。
2. 再用 `rows_need_siso_before_mux` 看这些峰值在经过 hybrid hard-finish 过滤后，还剩多少真正的 Chase/soft 需求。
3. 最后结合 `rows_unscheduled` 判断这些需求峰值是否真的引发了资源不足。

如果三者同时出现下面现象：

- `need_decode_after_early_stop` 有尖峰
- `rows_need_siso_before_mux` 也有对应尖峰
- `rows_unscheduled` 在这些时段显著增大

那么就可以比较有把握地说：

> 当前系统里确实存在“瞬时 soft/Chase 需求 burst”，而且它已经开始转化成资源不足症状。

这时再继续推进双虚拟流离线重放和共享资源池建模，就是有依据的。

### 3.1.5 这一阶段的通过标准

若观察到：

- `rows_need_siso_before_mux` 存在明显尖峰
- `rows_unscheduled` 主要集中在少数 invocation / tile
- BER burst 与这些尖峰时段有相关性

则说明“双流削峰”这个方向值得继续做。

如果看不到这些现象，则说明问题可能不是共享 Chase 峰值，而是别的因素，后续双流方案的优先级应降低。

## 3.2 第二阶段：基于单流 trace 做离线双流重放

### 3.2.1 目标

这一阶段回答：

> 即使不改主链路，只把当前单流时序样本离线重排成两路，是否已经能看到统计复用收益。

这是整个方案最关键的“成立性验证”。

### 3.2.2 核心思想

从第一阶段拿到单流时序样本后，不马上真的实现两路窗口，而是在离线分析里先构造：

- 虚拟流 A 的请求序列
- 虚拟流 B 的请求序列

然后观察这两路流在共享资源池视角下叠加后的总需求序列：

```text
a_total(t) = a_A(t) + a_B(t)
```

这里的关键不再是引入一个额外的可调“相位偏移量”，而是：

- 先构造两路相互独立的虚拟请求流
- 再观察这两路流在自然推进下叠加后的统计行为

这样更贴近当前可实现的硬件语义：  
两路流各自独立维护窗口和控制状态，只在共享 Chase 资源层面汇合，而不是要求实现一个额外可调的时间偏移量。

这里的 `a(t)` 在当前软件验证中，可以直接用：

```text
a(t) = rows_need_siso_before_mux
```

必要时也可以进一步定义：

```text
a_effective(t) = rows_need_siso_before_mux - rows_hard_finish_adjust
```

但第一版建议先不要复杂化，直接先用现有统计量。

### 3.2.3 如何构造两路虚拟流

建议从最简单的两种构造方式开始：

#### 方案 A：按 invocation 奇偶切分

- A 流取奇数 invocation
- B 流取偶数 invocation

优点：

- 实现简单
- 不需要理解更多内部结构

缺点：

- 和真实“独立窗口流”未必完全一致

#### 方案 B：按 tile/window 时序块切分

例如：

- 连续若干 invocation 作为 A
- 下一段作为 B

这种方式更接近“两个相互独立推进的流”，但实现略复杂。

第一轮建议先用方案 A，快速看有没有削峰潜力。

### 3.2.4 建议统计指标

对不同的双流构造方式，至少统计：

- `max(a_total)`
- `p99(a_total)`
- `p99.9(a_total)`
- `p99.99(a_total)`
- 高于给定 budget 的超限次数
- 连续超限长度分布

预算可以直接取当前 tile 的 `SISO_ACTIVE_LIST[t]`，或者先统一取一个标称预算值。

### 3.2.5 这一阶段的通过标准

如果某种双流构造方式下，观察到：

- 峰值明显下降
- 高分位数明显下降
- 超限持续长度明显缩短

则说明双流独立请求共享在统计上是有潜力的。

如果不同双流构造方式下都几乎没有改善，说明两路请求高度相关，这个方案就不值得继续做更重的工程接入。

## 3.3 第三阶段：软件级共享 Chase 资源池调度仿真

### 3.3.1 目标

第二阶段只回答“统计上会不会削峰”，第三阶段回答：

> 在真的存在共享资源池、队列和有限服务能力时，系统是否还能满足时序与质量目标。

这一阶段开始引入“有限资源”的概念，但仍然不要求改主解码流程。

### 3.3.2 建议建立的最小模型

建议先不要直接建复杂 RTL 模型，而是做一个离散事件 / 周期级软件调度器。

最小模型包含：

- 两路请求输入序列
- 一个共享 Chase 服务池，容量为 `C`
- 一个 FIFO 或优先队列
- 每个请求的服务时长 `service_cycles`
- 每个请求的最晚完成时刻 `deadline`

### 3.3.3 当前代码里如何定义请求

在当前软件框架下，一条“请求”不必一开始细到 row 级全部真实元��据，可以先用近似模型：

- `stream_id`
- `invocation`
- `tile_index`
- `demand = rows_need_siso_before_mux`

然后把一个 tile 的 demand 拆成多个单位请求，或者把它看成一个 batch 请求都可以。

为了贴近现有实现，我更建议第一版先做“按行展开”的单位请求模型：

- 一个需要 soft decode 的 row = 一个 Chase 请求

这样更容易和 `siso_budget`、`rows_unscheduled` 对齐。

### 3.3.4 deadline 如何定义

这是这一阶段最需要补的工程约束。

第一版建议不要引入复杂真实时钟，而是使用“相对 deadline”近似：

- 某个 invocation 生成的请求，必须在 `D` 个调度周期内完成
- 或必须在下一个 window / tile 写回点之前完成

这里的 `D` 可以先用参数化方式扫：

- `D = 1 tile`
- `D = 2 tiles`
- `D = 1 window step`

先看结论对 `D` 是否敏感。

### 3.3.5 调度策略建议

第一版建议分三层验证，不要一开始就上复杂优先级：

1. `FIFO`
2. `age-first`
3. `deadline-first`

若这三种简单策略都无法获得收益，再上 reliability-aware 调度才有意义。

### 3.3.6 建议输出指标

- 最大队列深度
- 平均队列深度
- deadline miss 次数
- deadline miss 比例
- 超限 burst 长度
- 同等 miss 约束下所需 Chase 容量 `C`

### 3.3.7 这一阶段的通过标准

如果在合理的 `C` 和 `D` 下：

- 队列可控
- deadline miss 很低
- 相比单流峰值配置，所需 Chase 容量下降明显

那这个方案就进入“值得做主链路软件接入验证”的阶段。

## 3.4 第四阶段：最小侵入式接入现有主链路

### 3.4.1 目标

前三阶段都成立后，再进入：

> 在当前 C++ 仿真主链路中，最小侵入地实现“双流独立状态 + 共享 Chase 调度”的软件验证版本。

这一步的目标不是做最终架构，而是验证：

- 两路独立窗口状态是否可在当前代码组织下跑通
- 共享 Chase 调度是否会改变 BER / burst 行为
- 当前代码结构是否适合承载这种调度层

### 3.4.2 不建议一开始就怎么做

不建议第一版就：

- 重写 `pipeline_runner`
- 改动现有 `ofec_decode_llr_*` 大量内部逻辑
- 把所有 tile / row 处理都并成统一多流框架

这样风险太大，而且很难定位问题。

### 3.4.3 建议的最小实现方式

建议新增一个单独的实验 runner，例如：

- `apps/ofec_dual_stream_probe.cpp`
- `src/ofec_dual_stream/...`

让它做下面几件事：

1. 从同一批输入中构造两路虚拟流状态。
2. 各自维护独立的 tile/window 进度。
3. 在“本该进入 Chase 的位置”不直接本地解码，而是先生成请求。
4. 请求统一进入一个软件共享池。
5. 共享池调度后，再把结果回写到各自流的状态里。

### 3.4.4 当前代码里最适合的挂钩位置

从当前代码结构看，最适合挂实验层的位置不是最底层 Chase 实现，而是 tile 调度层附近。

原因是：

- [ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp:331) 已经有“当前 tile 的 soft candidate 需求量”和“unscheduled 结果”
- [ofec_tile_dispatch.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_dispatch.ipp:1) 已经把 row 分成 `HardFinish / SoftDecode / EarlyStopAction / Unscheduled`

因此，一个比较现实的第一版接法是：

- 保持每路流在 early-stop / hybrid prepass 之前都独立运行
- 在 `run_mux_on_soft_candidates(...)` 之前或附近，截获两路的 soft candidate 集合
- 用一个共享调度器替代当前单 tile 本地 budget 裁剪

这样做的好处是：

- 不用一开始就去改 Chase 内核
- 先复用当前最清楚的“谁需要 Chase”语义
- 先验证“共享调度”是否真能减少 `rows_unscheduled`

### 3.4.5 第一版可以先不做的事情

为了控制复杂度，第一版建议暂时不做：

- 复杂 reliability-aware scheduler
- 多级 Chase 降级
- 完整 RTL 风格回写 crossbar 模型
- 任意多流扩展

只做：

- 两流
- 固定共享预算
- 简单调度策略

足够验证主想法了。

## 4. 每个阶段的交付物建议

为了让整个验证路径更可执行，建议每个阶段都产出明确交付物。

### 第一阶段交付物

- 单流 invocation 级时序 CSV
- `rows_need_siso_before_mux` 时序图
- `rows_unscheduled` 时序图
- BER burst 对齐分析图

### 第二阶段交付物

- 离线双流重放脚本
- 不同双流构造方式下的峰值 / 分位数对比表
- 双流重放前后需求分布对比图

### 第三阶段交付物

- 共享 Chase 资源池软件仿真器
- 容量 `C` 扫描结果
- queue depth / deadline miss 报表

### 第四阶段交付物

- 新增实验 runner
- 双流共享调度软件验证结果
- 和单流基线的 BER / burst / unscheduled 对比

## 5. 推荐的实施顺序

建议的真实执行顺序如下：

1. 先增强现有单流 trace 导出。
2. 做离线重放脚本，不改主链路。
3. 做共享池调度仿真器，不改主链路。
4. 只有前 3 步结果明显成立，才新增双流实验 runner。

这个顺序能最大限度复用现有代码，同时避免过早进入大改。

## 6. 总结

基于当前仓库，要验证双虚拟流共享 Chase 方案，不需要一开始就重写解码框架。  
更现实、风险更低的路径是：

1. 用现有 tile 统计先确认单流确实有瞬时 burst。
2. 用现有 trace 做离线双流重放，先验证是否能削峰。
3. 再加入共享资源池调度仿真，验证队列和 deadline 行为。
4. 最后才做最小侵入的双流软件接入。

如果前两阶段就看不到明显收益，那就说明这条路线不值得继续重投入；  
如果前三阶段都成立，再进入软件主链路接入就会更有把握，也更容易说服后续硬件实现。
