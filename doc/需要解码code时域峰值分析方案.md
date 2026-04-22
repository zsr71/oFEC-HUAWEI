# 需要解码 code 时域峰值分析方案

## 1. 问题背景

当前已经可以从 probe 数据里拿到每次 early-stop 检查后的统计量。

如果把“需要解码的 code”定义为：

```text
need_decode = 32 - rows_passed
```

那么在时域上通常会观察到这样一种现象：

- `need_decode` 在整个解码过程中不是平滑变化的
- 而是会出现多个峰值
- 其余大部分时间都比较低

你关心的是：

- 这种峰值意味着什么
- 它是否影响最终 BER
- 如何定量评估它的影响
- 是否需要进一步把“pre-BER / post-BER”也画成时域曲线

## 2. 现象的基本理解

这种峰值现象一般不能直接理解成“坏事”或者“好事”，它更像是：

- 解码状态在某些时间段集中进入困难区
- 某些 tile 或某些 window 的有效调度压力突然上升
- 译码器在局部时段需要处理更多未早停的 row

换句话说，`need_decode` 的峰值本身只说明：

- 当前时刻需要正常解码的 code 数量很多

它并不自动等于：

- 最终 BER 一定更差
- 或者系统一定不稳定

真正要判断的是：

- 峰值是否和错误率上升同步
- 峰值持续多长
- 峰值之后是否有恢复
- 峰值是否反复出现在相近的时间位置

## 3. 推荐的分析思路

建议把分析拆成三层。

### 3.1 第一层：只看 `need_decode` 的时域形态

最直接的图就是：

- 横轴：`invocation` 或时间序号
- 纵轴：`need_decode = 32 - rows_passed`

这张图的目标不是解释 BER，而是先找出：

- 峰值出现在哪些区间
- 峰值是否成簇
- 峰值是否有周期性
- 峰值持续多长

可以进一步做：

- 单次曲线
- 滑动平均曲线
- 峰值检测标记

### 3.2 第二层：把峰值和误码统计对齐

真正有价值的是把 `need_decode` 和 BER 放到同一条时间轴上看。

如果你的程序里已经有 per invocation 的样本，那么就可以同时画：

- `need_decode(invocation)`
- `pre_errs(invocation)` / `post_errs(invocation)`
- `pre_ber(invocation)` / `post_ber(invocation)`

这样可以看出：

- 峰值是否先于 BER 上升
- BER 是否在峰值期间更高
- 峰值消失后 BER 是否恢复

### 3.3 第三层：做分段统计，而不是只看单点曲线

建议把时间轴切成几个区间：

- 峰值前
- 峰值期间
- 峰值后

然后分别统计：

- 平均 `need_decode`
- 平均 `pre_ber`
- 平均 `post_ber`
- 平均 `unscheduled`
- 峰值持续长度

如果峰值期间 BER 明显更高，说明该区间确实在拉低性能。

如果峰值期间 BER 没有明显变化，那它可能只是内部调度压力变化，不一定是性能问题。

## 4. 如何评估影响

这里建议用三个指标。

### 4.1 相关性

看 `need_decode` 和 BER 是否相关。

例如：

- `need_decode` 高的时候，`post_ber` 是否也高
- `need_decode` 低的时候，`post_ber` 是否也低

这可以先做一个简单散点图：

- 横轴：`need_decode`
- 纵轴：`post_ber` 或 `post_errs`

如果散点图显示明显正相关，说明 `need_decode` 峰值和性能退化有关。

### 4.2 峰值持续时间

单个峰值不一定有影响。

更值得关注的是：

- 峰值是否持续 40~50 个 invocation
- 峰值期间是否覆盖了足够多的比特

持续时间越长，越可能对最终 BER 产生可见影响。

### 4.3 峰值对累计 BER 的边际贡献

可以分别计算：

- 只统计峰值前的累计 BER
- 只统计峰值期间的累计 BER
- 只统计峰值后的累计 BER

如果峰值期间的累计 BER 显著高于前后两段，那么它就是影响性能的重要区间。

## 5. 关于“画 pre/post BER 随时间变化”

你的想法是对的，这个图很有价值。

但要注意一个事实：

- `BER` 本身是需要一定统计量的
- 单个 invocation 的 BER 可能会很抖
- 尤其在低 BER 区间，单点经常是 0 error

因此建议不要只画单点 `BER`，而是画：

- 单点 BER
- 滑动平均 BER
- 累计 BER

更稳的做法是：

```text
rolling_BER(invocation, W)
```

也就是用一个固定窗口 `W` 做滑动平均。

这样可以同时看到：

- 短期峰值
- 以及峰值是否真的对应 BER 上升

## 6. 程序里怎么实现

### 6.1 最小实现

第一版可以只记录：

- `invocation`
- `tile_index`
- `rows_passed`
- `need_decode = 32 - rows_passed`

然后在 Matlab 里做：

- 画 `need_decode` 时域曲线
- 找峰值区间

这一步最容易落地。

### 6.2 进阶实现

如果要把 BER 也带上，建议再记录：

- `pre_errs`
- `pre_bits`
- `post_errs`
- `post_bits`

这样就能在同一条时间线上画：

- `need_decode`
- `pre_ber`
- `post_ber`

### 6.3 更完整的实现

如果你想更准确地分析“哪个峰值影响了最终 BER”，那最好记录每个 invocation 对应的：

- `tile_index`
- `window_idx`
- `rows_passed`
- `need_decode`
- `pre_errs`
- `post_errs`
- `pre_bits`
- `post_bits`

这样后面可以直接在程序里做分段统计，不必全靠 Matlab 手工切。

### 6.4 两种实现方案对比

这里其实有两条路可以走。

**方案 A：运行时直接记录 tile 级 BER 样本**

做法是：
- 每次 tile 处理完，就把该 tile 的误码数、比特数记下来
- 最后导出一份按 `invocation + tile_index` 排列的 CSV

优点：
- 和 `need_decode` 的时域样本天然同粒度
- 后面 Matlab 很容易直接做对齐
- 能保留每个 tile 的局部波动

缺点：
- 需要在 tile 执行链路里额外埋点
- 如果想记录很多中间量，代码会更散一点

如果真要落地，通常要改：
- [apps/ofec_ber_window_probe.cpp](/home/zsr71/projects/newcode/apps/ofec_ber_window_probe.cpp)
  - 增加 tile 级 BER 样本导出
- [include/newcode/pipeline_runner.hpp](/home/zsr71/projects/newcode/include/newcode/pipeline_runner.hpp)
  - 如果希望 probe 结果里直接带 tile 级 BER 样本，需要扩 `PipelineResult`
- [src/common/pipeline/pipeline_runner.cpp](/home/zsr71/projects/newcode/src/common/pipeline/pipeline_runner.cpp)
  - 在 pipeline 结束后把 tile 级统计挂到结果对象上
- [src/rx/ber/ber.cpp](/home/zsr71/projects/newcode/src/rx/ber/ber.cpp)
  - 如果要新增按 tile 切分的 BER 统计函数，可以在这里补

**方案 B：所有 window / tile 都跑完之后，再统一做 tile 级 BER 统计**

做法是：
- 运行时只保留 tile 的原始误差计数或必要的中间结果
- 程序结束后再统一汇总成 tile 级 BER

优点：
- 主译码链路更干净
- 对现有逻辑侵入更小
- 更适合先做第一版验证

缺点：
- 如果没有提前保存足够的原始信息，后面很难严格反推出 tile 级 BER
- tile 的时间对齐能力不如方案 A 直接

如果真要落地，建议新加一个 tile 版 BER 统计函数，而不是直接改旧函数：
- 保留现有 `compute_ber_per_window(...)`
- 新增一个例如 `compute_ber_per_tile_window(...)` 的函数
- 输入仍然是一维比特流，但窗口长度换成 tile 尺度
- 调用点可以放在新的 probe app 或 pipeline 结果后处理里

**我的判断**

如果你的目标是：
- 先分析最后一个 tile 的 `need_decode` 峰值和 BER 的相关性
- 并且希望后面 Matlab 里更容易对齐

那更推荐方案 A。

如果你的目标是：
- 先尽量少改代码
- 只要能得到 tile 级 BER 的汇总结果就行

那方案 B 更保守。

更实际的折中是：
- 先做方案 B，保证能算出 tile 级 BER
- 如果后面发现时序对齐不够，再补方案 A 的细粒度样本

## 7. 实现难度判断

你的感觉是对的：

- **只画 `need_decode` 时域图，不难**
- **把 BER 也按时间对齐，难度中等**
- **把 pre/post BER 严格做成按 invocation 的时域统计，改动会明显大一些**

为什么会这样：

- 现在的 BER 统计主要是“整次运行后汇总”
- 而你想看的，是“每个时间点的局部统计”
- 这意味着要把原来只在结束时输出的统计，拆到运行过程中逐步记录

所以第一版最现实的做法是：

1. 先把 `need_decode` 的时域图画出来
2. 再看是否需要把 `pre/post BER` 按 invocation 补出来

这比一开始就把所有 BER 时序接口都重构掉更稳。

## 8. 推荐结论

如果你的目标是分析“最后一个 tile 的 `need_decode` 在时间上多峰值变化的影响”，我建议按下面顺序来：

1. 先画 `need_decode = 32 - rows_passed` 的时域曲线
2. 再做滑动平均，确认峰值是否稳定存在
3. 再把峰值区间和 `post_ber` 对齐
4. 如果确实相关，再决定是否要在程序里加入 invocation 级的 `pre/post BER` 时序输出

一句话总结：

- `need_decode` 的峰值值得分析
- 最先要做的是确认它和 BER 是否相关
- 程序里可以实现，但第一版建议先做轻量统计，再逐步扩展到 BER 时序
