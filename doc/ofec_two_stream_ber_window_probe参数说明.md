# `ofec_two_stream_ber_window_probe` 参数说明

本文针对计划新增的 two-stream probe 入口，整理一份参数说明文档。

这里的目标不是直接讨论某一版具体代码怎么写，而是先把参数按类别拆清楚，说明：

- 这类参数是什么
- 这类参数控制的是哪一层行为
- 它对最终导出的 probe 数据会产生什么影响
- 哪些参数第一版必须有
- 哪些参数可以先保留但暂时不重点使用

为避免歧义，本文里的 “two-stream probe” 指的是一个未来新增的、风格接近
[apps/ofec_ber_window_probe.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_ber_window_probe.cpp)
的双流版探针入口，暂记为：

- `apps/ofec_two_stream_ber_window_probe.cpp`

它的定位是：

- 不做大范围 low-BER 扫描
- 也不只看单点总 BER
- 而是针对 two-stream shared 主流程，导出 A/B 两路在 window / tile / seed 维度上的细粒度统计

---

## 1. 参数设计的总原则

two-stream probe 的参数，建议按下面三条原则组织。

### 1.1 保持和 one-shot two-stream app 对齐

probe 最重要的作用之一，是帮助解释
[apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp)
里观察到的现象。

所以 probe 的大部分主流程参数，应该默认和 one-shot app 一致，包括：

- decoder 路径
- Chase 参数
- early-stop 参数
- hybrid 参数
- MUX 参数
- alpha / beta / SISO 预算

这样 probe 看到的时域分布，才能真正解释 one-shot app 的总 BER 行为。

### 1.2 把“前端独立参数”和“共享解码参数”分开

two-stream 的特殊点在于：

- A/B 两路前端是独立生成、独立过信道、独立量化的
- 但后面的 soft tile 路径会进入一个 shared 64-code 域

所以参数上一定要区分：

- 哪些参数作用在 A/B 各自前端
- 哪些参数作用在 merged shared 域

否则后面看 probe 结果时，会很难判断某个现象到底来自：

- 两路前端差异
- 还是 shared 调度差异

### 1.3 第一版先服务“定位问题”，不是服务“做全功能平台”

probe 和 sweep 不一样。

sweep 主要用于：

- 看 BER 曲线
- 做点级统计收敛

probe 主要用于：

- 看时域分布
- 看 seed 间复现性
- 看 A/B 是否对称
- 看 shared 域是不是在某些 tile 上偏流

因此第一版参数不需要做成“大而全的扫描平台”，而应该优先保证：

- 参数定义清晰
- 输出结构稳定
- 能快速复现并定位时域问题

---

## 2. 参数总览

建议把 two-stream probe 的参数分成下面几类：

1. 输出与运行组织参数
2. `Eb/N0` 与 seed 参数
3. A/B 前端生成参数
4. 信道与量化参数
5. 共享解码主参数
6. early-stop 参数
7. alpha / beta / SISO 预算参数
8. MUX 调度参数
9. Hybrid 参数
10. probe 专用观测参数
11. Debug / trace 参数

下面逐类展开。

---

## 3. 输出与运行组织参数

这一类参数不直接改变解码算法本身，而是决定：

- 输出文件放哪里
- 任务怎么并行
- probe 是按多少个 seed、多少个 `Eb/N0` 组织起来的

这一类参数很像
[apps/ofec_ber_window_probe.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_ber_window_probe.cpp)
顶部的输出与并行参数。

### 3.1 `kLabelPrefix`

定义：

- 运行标签前缀

作用：

- 进入日志文件名
- 进入内部每个任务的 label
- 便于区分不同 probe 入口或不同实验方案

建议：

- 第一版必须保留
- 应明确体现 two-stream 身份，例如 `two_stream_ber_window_probe`

### 3.2 `kOutputDir`

定义：

- probe 所有输出文件的根目录

作用：

- 控制 CSV、日志文件的落盘位置

建议：

- 第一版必须保留
- 最好单独目录，例如 `data/two_stream_ber_window_probe`

### 3.3 `kQuietPipeline`

定义：

- 是否压低单次 pipeline 内部的控制台输出

作用：

- `true` 时能避免 probe 跑多个 seed 时日志爆炸
- `false` 时更适合调试单个异常 seed

建议：

- 第一版保留
- 默认建议 `true`

### 3.4 `kMaxParallelEbN0`

定义：

- 外层 `Eb/N0` 任务的最大并行度

作用：

- 控制是否同时并行跑多个 `Eb/N0`

建议：

- 第一版保留
- 如果第一版只打算固定一个 `Eb/N0` 点，这个参数仍然可以保留，但不一定常用

### 3.5 `kMaxParallelSeeds`

定义：

- 单个 `Eb/N0` 点内部 seed 任务的最大并行度

作用：

- 控制同一 `Eb/N0` 下多个 seed 是否并行执行

建议：

- 第一版必须保留
- 因为 probe 的一个重点就是看 seed 间稳定性

---

## 4. `Eb/N0` 与 seed 参数

这一类参数决定“要在哪些条件下采样”。

probe 不像 sweep 那样追求大范围扫描，但依然需要：

- 指定测试点
- 指定重复 seed 数量

### 4.1 `kEbN0List`

定义：

- 要测试的 `Eb/N0` 列表

作用：

- 决定 probe 在哪些噪声强度下采样

建议：

- 第一版建议只放少量点
- 通常 1 到 3 个点就够了

原因：

- probe 更偏定位工具，不需要像 sweep 那样扫很多点
- 过多点会稀释分析精力

### 4.2 `kBitgenSeedBaseA` / `kChannelSeedBaseA`

定义：

- Stream A 的基础信息比特 seed / 基础信道 seed

作用：

- 控制 A 路每个 seed 任务的随机源

### 4.3 `kBitgenSeedBaseB` / `kChannelSeedBaseB`

定义：

- Stream B 的基础信息比特 seed / 基础信道 seed

作用：

- 控制 B 路每个 seed 任务的随机源

建议：

- two-stream probe 不应该继续沿用单流 probe 只有一套 base seed 的写法
- 应显式分成 A/B 两组

原因：

- 两路原生独立是 two-stream 的核心前提
- 如果参数层面不区分 A/B 基础 seed，后面很难做 A/B 交换实验

### 4.4 `kSeedCount`

定义：

- 每个 `Eb/N0` 点要跑多少组 seed

作用：

- 控制 probe 结果在 seed 维度上的重复次数

建议：

- 第一版必须保留
- 如果只是结构对通，可以先用很小值
- 如果要看时域偏差是否稳定，建议至少大于 1

---

## 5. A/B 前端生成参数

这一类参数作用在每一路各自的发端生成阶段。

它们本身不涉及 shared 域，但会影响：

- 每一路进入 shared path 之前的原始统计
- A/B 前端是否天然有差异

### 5.1 `kGenerateRandomBits`

定义：

- 是否发送随机信息比特

作用：

- `true` 时每个 seed 发送随机信息序列
- `false` 时发送全 0 信息

建议：

- 第一版保留
- 默认建议 `true`

原因：

- two-stream 的真实行为更接近随机比特
- 全 0 更适合极简调试，不适合长期作为 probe 主设置

### 5.2 `kInterleaverName`

定义：

- 交织器名称

作用：

- 控制每一路发端/收端使用的交织配置

建议：

- 第一版必须与 one-shot app 对齐
- 通常先保持 `identity`

### 5.3 `kBitsPerSymbol`

定义：

- 每个调制符号携带的比特数

作用：

- 决定 BPSK / QAM 等调制入口

建议：

- 第一版保留
- 默认与当前 two-stream one-shot app 一致

---

## 6. 信道与量化参数

这一类参数决定两路前端过信道后的 LLR 口径，以及进入 shared 解码前的量化表现。

### 6.1 `kLlrBits`

定义：

- LLR 位宽

作用：

- 控制是浮点直通还是量化表示

影响：

- pre-FEC quantized hard BER
- 量化饱和率
- 进入 shared 域的实际 LLR 离散程度

建议：

- 第一版必须和 one-shot two-stream app 对齐

### 6.2 `kQuantClipRatio`

定义：

- 动态 clip 比例

作用：

- 控制量化前 clip 的计算口径

two-stream 特别说明：

- 这里最终影响的不是单路 clip，而是 shared 量化口径
- 当前 two-stream shared 实现里，会基于两路前端共同算出一个 shared quant clip

建议：

- 第一版必须保留
- 因为它直接影响 A/B 两路的量化一致性

---

## 7. 共享解码主参数

这一类参数决定 shared decoder 走哪条主路径。

它们是 two-stream probe 和单流 probe 最大的区别之一。

### 7.1 `kDecoderName`

定义：

- 当前 probe 使用的 shared decoder 路径名称

作用：

- 决定 two-stream shared runner 底层选用哪种解码核心

建议：

- 第一版必须保留
- 默认应与 one-shot app 当前配置一致

### 7.2 `kChaseL`

定义：

- Chase L 参数

作用：

- 决定最不可靠位展开规模

### 7.3 `kChaseNTestOverride`

定义：

- Chase 测试 pattern 数覆盖值

作用：

- 若为负，则通常按 `2^L`
- 若为非负，则强制使用指定测试数

### 7.4 `kChaseTopkKeep`

定义：

- Top-K / pruned 路径保留的候选数

作用：

- 影响候选保留规模

### 7.5 `kChaseGroupMinimaBits`

定义：

- group minima 路径的分组 bit 数

作用：

- 决定分组最小值策略的粒度

### 7.6 `kNormalizeExtrinsic`

定义：

- 是否对 decoder 输出的 extrinsic 做归一化

作用：

- 影响 tile 输出回写后的外信息幅度口径

### 7.7 `kNormalizeKnownPrefixTail`

定义：

- 是否对 known-prefix 后的尾部 LLR 做归一化

作用：

- 影响前后端边界处的 LLR 处理口径

建议：

- 以上这一组参数，第一版都应该完整保留
- 即使某些 decoder 路径暂时没有真正用到，也建议保留参数区完整性

---

## 8. early-stop 参数

这一类参数控制 shared 64-code 域上 soft tile 的 early-stop 判定与动作。

它们对 probe 的意义非常大，因为：

- 它们直接影响哪些 row 不再进入后续 soft competition
- 也会影响某些 tile 上 A/B 是否开始偏流

### 8.1 `kEnableEarlyStop`

定义：

- early-stop 总开关

作用：

- `true` 时启用 early-stop
- `false` 时所有 row 都继续进入后续路径

### 8.2 `kEarlyStopEnableList`

定义：

- 按 tile 覆盖 early-stop 开关的列表

作用：

- 控制每个 tile 是否启用 early-stop

### 8.3 `kEarlyStopConditionMode`

定义：

- early-stop 条件模式

作用：

- 决定使用哪类条件

### 8.4 `kEarlyStopConditionModeList`

定义：

- 按 tile 覆盖条件模式的列表

作用：

- 允许不同 tile 使用不同判定逻辑

### 8.5 `kEarlyStopActionMode`

定义：

- 命中 early-stop 后的动作模式

作用：

- 决定命中后是如何生成输出 LLR / 外信息

### 8.6 `kEarlyStopActionModeList`

定义：

- 按 tile 覆盖动作模式的列表

作用：

- 允许不同 tile 使用不同动作策略

### 8.7 `kEarlyStopBindGroupSize`

定义：

- 条件模式 1 下的组绑定大小

作用：

- `1` 表示逐 row 判定
- `4` 表示按组绑定通过

### 8.8 `kEarlyStopBindGroupSizeList`

定义：

- 按 tile 覆盖绑定组大小的列表

### 8.9 `kEarlyStopCondV1RequireBch`

定义：

- v1 条件下是否要求 BCH syndrome 通过

### 8.10 `kEarlyStopCondV1RequireOverall`

定义：

- v1 条件下是否要求 overall parity 通过

### 8.11 `kEarlyStopV2LlrAbsThreshold`

定义：

- v2 条件下定义“不可靠位”的 `|LLR|` 阈值

### 8.12 `kEarlyStopV2MaxUnreliableBits`

定义：

- v2 条件允许的不可靠 bit 数上限

### 8.13 `kEarlyStopCondV2IncludeOverall`

定义：

- v2 条件下是否把 overall parity bit 计入不可靠位统计

### 8.14 `kEarlyStopActionResidualDivisor`

定义：

- 某些动作模式中 residual 的缩放除数

### 8.15 `kEarlyStopActionHardLlrMag`

定义：

- 硬完成动作输出的固定 `|LLR|` 幅度

建议：

- 第一版必须完整保留这组参数
- 因为 probe 的一个重要用途，就是看 early-stop 是否在某些 tile 或某一路上出现系统性偏差

---

## 9. alpha / beta / SISO 预算参数

这一类参数控制共享域里每个 tile 的软信息幅度和 soft competition 预算。

在 two-stream 场景里，它们通常是最敏感的一组参数之一。

### 9.1 `kAlphaExplicit`

定义：

- 每个 tile 的 alpha 显式列表

作用：

- 控制某些软信息处理环节里的缩放系数

### 9.2 `kBetaExplicit`

定义：

- 每个 tile 的 beta 显式列表

作用：

- 控制每个 tile 的 Chase / fallback 输出幅度口径

### 9.3 `kEarlyStopActionBetaExplicit`

定义：

- 每个 tile 专门用于 early-stop 动作的 beta 列表

作用：

- 允许 early-stop 输出和普通 soft path 使用不同幅度

### 9.4 `kSisoActiveList`

定义：

- 每个 tile 在 shared 域里的 SISO 预算列表

作用：

- 直接控制每个 tile 能有多少 merged rows 进入 soft 竞争或 soft 核

two-stream 特别说明：

- 在单流里它是单路 row 预算
- 在 two-stream shared 里，它作用在 merged 的 64-code 域上

这意味着：

- `32` 不再表示“一路 32 个”
- 而是“共享域总共可用 32 个”

这是 two-stream probe 非常需要观测的一点。

建议：

- 第一版必须保留
- 而且后续 shared probe CSV 里最好把它对应到实际 produced / failed 行为上

---

## 10. MUX 调度参数

这一类参数控制 shared 预算裁剪和 soft candidate 调度方式。

它们直接关系到：

- A/B 在共享域中如何竞争资源
- 是否存在顺序敏感
- 某些 tile 上会不会对某一路更偏置

### 10.1 `kMuxGroupG`

定义：

- MUX 分组数

作用：

- `1` 时通常表示全局共享预算池
- 大于 `1` 时表示把共享域切成若干组分别调度

### 10.2 `kMuxSchedulingMode`

定义：

- MUX 调度模式

作用：

- 控制 soft candidate 的裁剪规则

### 10.3 `kMuxPriorityRule`

定义：

- MUX 的优先级规则

作用：

- 在某些调度模式下，决定 row 的保留优先级

### 10.4 `kMuxEnableReconfig`

定义：

- 是否启用 reconfig 型 MUX 调度

作用：

- 决定是否走 staged / reconfig 调度逻辑

### 10.5 `kMuxExtraBypassEdges`

定义：

- reconfig 模式使用的旁路边集合

作用：

- 控制 MUX 图上的额外旁路连接

建议：

- 第一版应该完整保留
- 即使 probe 初版不导出所有 MUX 内部细节，至少参数必须完整记录

---

## 11. Hybrid 参数

这一类参数控制在 shared soft core 前，哪些 row 可以先走硬分类、硬完成或回填逻辑。

这部分也是 two-stream probe 很值得观察的一层，因为它常常会改变：

- 哪些 row 还会进入后续 shared soft competition
- A/B 两路在后续阶段的竞争结构

### 11.1 `kHybridEnable`

定义：

- hybrid 总开关

作用：

- 是��启用 hybrid 分流

### 11.2 `kHybridEnableList`

定义：

- 按 tile 覆盖 hybrid 开关

作用：

- 允许只在部分 tile 启用 hybrid

### 11.3 `kHybridHardLlrMag`

定义：

- hybrid hard-finish 输出的固定 `|LLR|` 幅度

### 11.4 `kHybridHardLlrMagList`

定义：

- 按 tile 覆盖 hard-finish 幅度的列表

### 11.5 `kHybridClassifierMode`

定义：

- hybrid 分类器模式

作用：

- 决定 row 如何被归类为：
  - 直接硬完成
  - 继续 soft
  - 或其它中间状态

### 11.6 `kHybridSisoBackfillMode`

定义：

- hybrid 的 SISO 回填策略

作用：

- 决定哪些分类结果还能重新进入 soft 预算竞争

### 11.7 `kHybridNormalizeSoftOnly`

定义：

- 是否只归一化 soft rows

作用：

- 控制 hybrid 后不同类型 row 的 LLR 口径一致性

建议：

- 第一版必须保留这一整组参数
- 因为你当前 two-stream one-shot app 已经明显依赖 hybrid 路径

---

## 12. probe 专用观测参数

这一类参数不改变 shared 解码逻辑本身，而是决定 probe 要导出哪些观察数据。

这部分是 two-stream probe 相比 one-shot app 最该单独设计的一块。

建议把 probe 专用观测参数再分成两层。

### 12.1 A/B 路时域 BER 导出开关

建议参数：

- `kDumpPerSeedPerWindowA`
- `kDumpPerSeedPerWindowB`
- `kDumpPerSeedPerTileA`
- `kDumpPerSeedPerTileB`
- `kDumpAggregatedWindowA`
- `kDumpAggregatedWindowB`
- `kDumpAggregatedTileA`
- `kDumpAggregatedTileB`

定义：

- 控制是否导出 A/B 各自的逐 seed、逐 window / tile 结果

作用：

- 服务于 burst 位置、window BER 分布、tile BER 分布分析

建议：

- 第一版必须有

### 12.2 shared 域样本导出开关

建议参数：

- `kDumpSharedTileSamples`
- `kDumpSharedCoreSummaryPerSeed`

定义：

- 控制是否导出 shared 域的 tile 级样本

作用：

- 服务于 A/B 竞争 shared 预算时的偏流定位

建议：

- 第一版可以先只保留参数定义
- 第二版再把 shared tile 样本 CSV 做完整

### 12.3 观测粒度开关

建议参数：

- `kRecordWindowStats`
- `kRecordTileStats`
- `kRecordEarlyStopSamples`
- `kRecordMuxSamples`
- `kRecordHybridSamples`

定义：

- 控制 probe 是否记录某类中间观察量

作用：

- 避免默认导出过多数据
- 便于在“轻量模式”和“深入排查模式”之间切换

建议：

- 第一版可以保留，但只实现其中最关键的几项

---

## 13. Debug / trace 参数

这一类参数主要服务于异常定位，不应该作为日常 probe 的主输出手段。

### 13.1 `kDecoderTraceEnable`

定义：

- decoder trace 总开关

作用：

- 是否把 trace 配置灌进 two-stream shared 主流程

### 13.2 `kDecoderTraceRow`

定义：

- 追踪目标 row

### 13.3 `kDecoderTraceCol`

定义：

- 追踪目标 col

### 13.4 `kDecoderTraceLogRead`

定义：

- 是否记录 tile 读取映射

### 13.5 `kDecoderTraceLogWrite`

定义：

- 是否记录 tile 写回映射

### 13.6 `kDecoderTraceLogMismatch`

定义：

- 是否记录写回不一致告警

建议：

- 第一版保留
- 默认关闭

原因：

- 这些参数更多是“显微镜”
- probe 本身应该先解决“望远镜”层面的问题，也就是时域统计结构

---

## 14. 哪些参数第一版必须实现

如果按“第一版先做 A/B 时域 BER probe”的思路，建议第一版必须实现下面这些类别：

1. 输出与运行组织参数
2. `Eb/N0` 与 seed 参数
3. A/B 前端生成参数
4. 信道与量化参数
5. 共享解码主参数
6. early-stop 参数
7. alpha / beta / SISO 预算参数
8. MUX 调度参数
9. Hybrid 参数

这些参数决定的是：

- probe 和 one-shot app 是否同口径
- probe 结果能不能真正解释当前 two-stream shared 行为

---

## 15. 哪些参数第一版可以只定义、不一定立即导出完整效果

下面这些参数，第一版可以先把参数位留好，但不一定一开始就把对应观测 CSV 做满：

1. shared 域样本导出开关
2. MUX 样本导出开关
3. Hybrid 样本导出开关
4. 深度 trace 参数

原因是：

- 第一版最需要先确认 A/B 的 window / tile BER 分布
- shared 内部样本可以在第二阶段逐步加

---

## 16. 参数之间最重要的几组关系

最后补一节，专门强调 two-stream probe 里最容易混淆的几组关系。

### 16.1 `A/B seed` 和 `shared 域偏流` 不是一回事

如果 A/B 两路基础 seed 不同，那么：

- A/B 前端本来就不完全相同

这时 probe 的作用是进一步区分：

- 差异是来自前端随机性
- 还是 shared 域把差异放大了

所以 seed 参数和 shared 参数必须分开记录。

### 16.2 `kSisoActiveList` 在 two-stream 中是 shared 预算

这是最关键的一点之一。

在 two-stream shared 里，`kSisoActiveList` 控制的是 merged 64-code 域的预算，而不是 A/B 各自独享的预算。

这意味着：

- 它天然可能导致 A/B 竞争
- 也天然需要 probe 去观察实际分配结果

### 16.3 early-stop / hybrid / MUX 是串联影响，不应孤立看

这三组参数不是互相独立的。

更准确地说：

- early-stop 先决定谁提前退出
- hybrid 再决定谁直接硬完成、谁还能参与竞争
- MUX 再对剩余 soft candidates 做预算裁剪

因此 probe 分析时不能只看其中一层参数，而要把三层放在同一个 tile 上一起看。

---

## 17. 建议的文档化方式

如果后面真的开始写
`apps/ofec_two_stream_ber_window_probe.cpp`，
我建议参数区就按本文的分组直接排：

1. 输出与运行组织参数
2. `Eb/N0` 与 seed 参数
3. A/B 前端生成参数
4. 信道与量化参数
5. 共享解码主参数
6. early-stop 参数
7. alpha / beta / SISO 预算参数
8. MUX 调度参数
9. Hybrid 参数
10. probe 专用观测参数
11. Debug / trace 参数

这样后面你自己看代码时，会比“按历史演化顺序堆参数”清晰很多。

---

## 18. 总结

two-stream probe 的参数设计，核心不是“再抄一份单流 probe 参数”，而是要把下面三层明确拆开：

- A/B 各自前端的独立参数
- 进入 merged 64-code 域之后的共享解码参数
- probe 自己为了导出时域统计而新增的观测参数

其中最关键的几组参数是：

- A/B 基础 seed
- `kLlrBits` / `kQuantClipRatio`
- `kDecoderName` / Chase 参数
- early-stop 参数
- `kSisoActiveList`
- MUX 参数
- Hybrid 参数

如果这几类参数定义清楚了，后面不管是做第一版的 A/B 时域 BER probe，还是第二版的 shared 域偏流 probe，结构都会比较稳。
