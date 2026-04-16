# `ofec_sweep3` 低 BER 扫描方案

## 1. 目标

当前 [ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 主要适合扫描常规 BER 曲线：

- 每个 `Eb/N0` 点对应一个独立 scenario
- 每个 scenario 由一个 worker 独立执行一次完整 Monte Carlo
- 最终直接读取这次运行的 `pre/post BER`

这种方式在纠后 BER 还比较高时足够直接，但当目标下降到 `1e-6`、`1e-7`、`1e-8` 甚至更低时，会遇到两个明显问题：

- 单个 `Eb/N0` 点需要极大量比特才能观察到足够多的 post-FEC errors
- 当前并行粒度是“按点并行”，高 BER 低点和低 BER 高点会一起拖慢整体扫描

因此建议**不改现有 `ofec_sweep` 逻辑**，而是新增一个专门面向低 BER 的 [ofec_sweep3.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep3.cpp)。

## 2. 为什么单独新增 `ofec_sweep3`

推荐单独做 `ofec_sweep3`，而不是在 `ofec_sweep` 上直接叠加新逻辑，原因是：

- `ofec_sweep` 现在的职责很清晰，就是“标准 scenario 扩展 + 每点一次完整仿真”
- 低 BER 扫描需要完全不同的运行策略：
  - 单点多 chunk 聚合
  - 自适应停止
  - 结果上界估计
  - 更细的 worker 调度
- 如果直接塞进 `ofec_sweep`，会让当前简单模式变复杂，后续维护和对比都不方便

因此更合适的结构是：

- `ofec_sweep`
  - 保持现有 BER 曲线扫描逻辑不变
- `ofec_sweep3`
  - 专门负责低 BER 区间扫描

这样后面更容易做 A/B 对比：

- 同一组参数
- `sweep` 看常规区间
- `sweep3` 看低 BER 区间

## 3. `ofec_sweep3` 的核心思路

### 3.1 基本思想

对于每个 `scenario + Eb/N0` 点，不再只跑一次大 Monte Carlo，而是拆成很多个较小的独立 chunk：

- 每个 chunk 使用一组独立 seed
- 每个 chunk 跑固定数量的信息比特
- 多个 chunk 可以并行跑
- 主线程持续聚合 chunk 的统计结果

最终每个点的 BER 不再来自“一次运行”，而是来自：

- `post_fec_errors_total = sum(chunk.post_fec.errors)`
- `post_fec_bits_total = sum(chunk.post_fec.total)`
- `ber = post_fec_errors_total / post_fec_bits_total`

### 3.2 为什么这比当前模式更适合低 BER

这样做有两个直接优势：

- 单个 `Eb/N0` 点可以同时吃掉多个 CPU worker，不再局限于一个点一个核
- 可以在达到统计目标后立即停止，而不是预先猜一个很大的 `NUM_INFO_BITS`

也就是说，`ofec_sweep3` 的并行粒度要从：

- 当前的“按 `Eb/N0` 点并行”

改成：

- “按 `scenario + Eb/N0 + chunk` 并行”

## 4. 推荐的停止条件

低 BER 扫描最重要的不是“每次跑多大”，而是“什么时候可以停”。

推荐第一版就支持下面三条停止条件。

### 4.1 目标错误数停止

当某个点累计到足够多的 post-FEC errors 后停止，例如：

- `target_post_errors = 50`
- 或 `100`

这是最直接、最稳的停止条件。

优点是：

- BER 估计的统计波动可控
- 容易和不同 `Eb/N0` 点公平比较

### 4.2 最大总比特数停止

为了防止某些极低 BER 点无限跑下去，需要设一个上限：

- `max_post_fec_total_bits`

例如：

- `1e8`
- `5e8`
- `1e9`

一旦累计比较的 post-FEC 总比特数超过这个值，就停止。

### 4.3 零错误上界停止

在极低 BER 区域，常见情况是：

- 已经跑了很多 bits
- 但 `post_fec.errors = 0`

这时不能简单报：

- `BER = 0`

更合理的是给一个上界，例如 95% 置信下常用近似：

- `upper_bound ~= 3 / total_bits`

因此可以增加一个停止准则：

- 若当前 `errors = 0`
- 且 `upper_bound` 已经低于目标 BER
- 则停止该点，并以“上界结果”记入输出

## 5. 推荐的运行模式

### 5.1 Chunk 模式

每个 chunk 用固定的信息比特数，例如：

- `chunk_num_info_bits = 2e6`
- 或 `5e6`

这样单个 chunk：

- 运行时间不会太长
- 结果粒度足够细
- 容易做动态调度

### 5.2 单点多 worker 聚合

在低 BER 扫描中，推荐对每个 `Eb/N0` 点采用：

- 一个聚合器
- 多个 worker chunk 并发

而不是继续用“每点一个 worker”的方式。

推荐支持：

- `max_parallel_chunks_per_point`

例如：

- `4`
- `8`
- `16`

### 5.3 高 BER 和低 BER 分开策略

后续如果要做得更完整，可以考虑双模式：

- 高 BER 区：
  - 沿用当前 `ofec_sweep` 风格
- 低 BER 区：
  - 使用 `ofec_sweep3` 的聚合扫描

但第一版不必把两者合并到一个 app 里，分别保留更清楚。

## 6. 结果输出建议

`ofec_sweep3` 不应只输出最终 BER，还应把低 BER 扫描的统计过程写出来，方便后处理。

推荐每个点输出：

- `scenario`
- `Eb/N0`
- `post_fec_errors_total`
- `post_fec_total_bits`
- `post_fec_ber`
- `ber_is_upper_bound`
- `ber_upper_bound`
- `chunks_completed`
- `stop_reason`
- `elapsed_seconds`

其中：

- `ber_is_upper_bound = 0/1`
- `stop_reason` 可选：
  - `target_errors`
  - `max_bits`
  - `zero_error_upper_bound`

这样后面你在 Matlab 或 Python 里画图时，不仅能画 BER 曲线，还能知道每个点的统计可信度。

## 7. 与当前代码库的对接方式

### 7.1 尽量复用现有的 `run_pipeline`

当前 [pipeline_runner.cpp](/home/zsr71/projects/newcode/src/common/pipeline/pipeline_runner.cpp) 已经能返回：

- `result.post_fec.errors`
- `result.post_fec.total`

这正好适合做 chunk 聚合。

所以 `ofec_sweep3` 第一版不需要改 decoder 主链，只需要：

- 每次构造一个小 `NUM_INFO_BITS` 的 `Params`
- 多次调用 `run_pipeline(...)`
- 把返回结果累加

### 7.2 尽量复用现有的 scenario 展开

当前 [ofec_sweep_scenarios.cpp](/home/zsr71/projects/newcode/src/ofec_sweep/ofec_sweep_scenarios.cpp) 已经能完成：

- 参数展开
- scenario 命名
- alpha/beta/early-stop/MUX 接线

因此 `ofec_sweep3` 最好不要重新发明一套配置体系，而是：

- 继续复用 `SweepParameterConfig`
- 继续复用 `build_scenarios(...)`
- 只改执行策略

这点和 `ofec_sweep2` 的思路类似：

- 配置体系复用标准 `sweep`
- 运行方式单独实现

### 7.3 新增聚合执行器

需要新增一个专门的执行器，例如概念上：

- `run_scenarios_low_ber_parallel(...)`

它和当前的 `run_scenarios_parallel(...)` 最大区别是：

- 当前版本：每个 scenario 只跑一次
- 新版本：每个 scenario 被拆成多个 chunk，持续聚合直到停止条件满足

## 8. 推荐新增的配置项

建议在 `ofec_sweep3` 顶层新增这些参数：

- `kChunkNumInfoBits`
- `kTargetPostErrors`
- `kMaxPostFecTotalBits`
- `kMaxParallelChunksPerPoint`
- `kEnableZeroErrorUpperBound`
- `kTargetBerUpperBound`
- `kConfidenceLevel`

这些参数只放在 `ofec_sweep3.cpp` 对应的配置区，不回灌到 `ofec_sweep.cpp`。

## 9. 推荐的实现顺序

### 9.1 第一步：先做最小可用版本

先实现：

- 复用 `build_scenarios(...)`
- 每个点拆成多个 chunk
- 聚合 `post_fec.errors / total`
- 停止条件只支持：
  - `target_post_errors`
  - `max_post_fec_total_bits`

这一版就已经能明显提升低 BER 扫描效率。

### 9.2 第二步：补充零错误上界

在第一版稳定后，再加：

- `errors = 0` 时的 BER 上界估计
- 对应 CSV 输出字段

这一步可以避免在超低 BER 点上跑得过久。

### 9.3 第三步：补充更细的调度和日志

后续再考虑：

- worker 资源在不同 `Eb/N0` 点之间动态分配
- chunk 进度日志
- per-point ETA
- 更丰富的 CSV

## 10. 总结

推荐方案不是修改现有 [ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)，而是新增一个专门面向低 BER 的 [ofec_sweep3.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep3.cpp)。

它的关键思想是：

- 不再把一个 `Eb/N0` 点交给一个 CPU 跑到底
- 而是把同一个点拆成多个独立 chunk 并行
- 持续聚合统计量
- 达到停止条件就停止

这样更适合扫到：

- `1e-7`
- `1e-8`
- 甚至更低

并且对现有代码库最友好，因为：

- decoder 主链基本不用动
- `SweepParameterConfig` 和 `build_scenarios(...)` 可以继续复用
- 风险主要集中在新的执行器层，而不是算法层
