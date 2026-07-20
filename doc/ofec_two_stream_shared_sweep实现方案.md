# `ofec_two_stream_shared_sweep.cpp` 实现方案

## 1. 目标

当前仓库里已经有两个很重要的参考入口：

- [apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp)
  - 负责跑单个固定参数点的 two-stream shared 主流程验证
- [apps/ofec_sweep3.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_sweep3.cpp)
  - 负责单流 low-BER 扫描，核心特点是“单点多 chunk 聚合 + 全局 worker 调度”

接下来需要的新入口 `ofec_two_stream_shared_sweep.cpp`，目标不是再发明一套新的 two-stream 解码主流程，而是：

- 复用当前已经接通的 two-stream shared runner
- 复用 `sweep3` 已经验证过的低 BER 扫描外壳
- 形成一个专门面向“两路原生独立流 + 共享 Chase 主流程”的扫描入口

这个 app 的定位可以概括为：

> 用 `sweep3` 的扫描外壳，去驱动 `two_stream_shared_runner` 的双流共享解码主流程。

第一版的重点不是一次把所有扫描维度都做全，而是先做一个：

- 只扫 `Eb/N0`
- 其它参数与 `ofec_two_stream_shared_decoder.cpp` 保持一致
- 能稳定输出 A/B 两路 BER 曲线的最小可用版本

---

## 2. 为什么要单独新增一个 app

不建议直接把 two-stream 逻辑塞进现有 [apps/ofec_sweep3.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_sweep3.cpp)，原因有三个。

### 2.1 单流 sweep3 和双流 shared sweep 的结果类型已经不同

`ofec_sweep3` 当前每个点只需要聚合一套：

- `pre_fec`
- `pre_fec_quantized_hard`
- `post_fec`

但 two-stream shared sweep 至少要同时输出：

- Stream A 的 `pre/post BER`
- Stream B 的 `pre/post BER`

也就是说，它已经不是“同一个结果结构多几列”这么简单，而是执行对象都变了。

### 2.2 two-stream shared sweep 的停止条件要按双流口径定义

单流时，一个点是否停止，只需要看这一条流的 post-FEC 累积统计。

双流时，如果只看其中一路，就会出现问题：

- A 已经累计到足够错误数
- B 还没有
- 但整个点被提前停掉

这样 B 的统计不够稳，后面的曲线会失真。

因此 two-stream shared sweep 需要一套自己的停止定义，而不是直接照搬单流版本。

### 2.3 two-stream shared sweep 的主目标仍然是 BER 曲线

这个 app 的主要用途，是把 two-stream shared 主流程放进 low-BER 扫描框架里，然后得到 A/B 两路的 BER 曲线。

像下面这些 shared 可观测性：

- A/B 是否公平竞争 shared 预算
- shared core 的 produced/failed 是否偏流
- quant clip 与饱和比例是否稳定

在单点入口
[apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp)
里依然很有价值，因为那边更适合逐点看内部行为。

但在 `sweep` 入口里，这些不是第一优先级。

也就是说，`ofec_two_stream_shared_sweep.cpp` 第一版更关心的是：

- shared 主流程在不同 `Eb/N0` 下的 BER 曲线是否合理
- A/B 两路曲线是否一致或符合预期
- sweep 结果与单点 app 能不能对上

因此这里直接定一个明确边界：

- `ofec_two_stream_shared_sweep.cpp` 第一版完全不聚合 shared 可观测性
- sweep 结果只服务于 BER 曲线分析
- 如果后面要看 shared 内部行为，统一回到单点 two-stream app 去看

因此更合适的做法仍然是：

- 保留 [apps/ofec_sweep3.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_sweep3.cpp) 作为单流 low-BER 扫描入口
- 新增 [apps/ofec_two_stream_shared_sweep.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_sweep.cpp) 作为双流共享扫描入口

这样职责会更清楚。

---

## 3. 与现有代码的关系

## 3.1 直接复用的部分

下面这些能力可以直接复用，不建议重写：

- [src/two_stream_shared/two_stream_shared_runner.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/two_stream_shared/two_stream_shared_runner.cpp)
  - 已经接通 two-stream shared 主流程
- [include/newcode/two_stream_shared_runner.hpp](/home/zsr71/projects/newcode_two_stream_shared_chase/include/newcode/two_stream_shared_runner.hpp)
  - 已经定义了 `Config` 和 `Result`
- [src/ofec_sweep/ofec_sweep_detail.hpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/ofec_sweep/ofec_sweep_detail.hpp)
  - 提供 `DualOut`、`build_ebn0_values(...)`、`resolve_worker_count(...)`
- [src/ofec_sweep/ofec_sweep_scenarios.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/src/ofec_sweep/ofec_sweep_scenarios.cpp)
  - 提供单流 sweep 场景展开逻辑

## 3.2 不建议直接复用的部分

下面这些不建议直接硬套：

- `ofec_sweep::run_sweep(config)`
- `ofec_sweep3` 现成的 `AggregatedPointResult`
- `ofec_sweep3` 现成的 `run_chunk(...)`

原因不是它们写得不好，而是它们默认的执行对象仍然是单流 `run_pipeline(...)`，结果类型和停止条件都按单流口径定死了。

因此更合适的方式是：

- 复用 `sweep3` 的“结构”
- 但自己在 `two_stream_shared_sweep.cpp` 里定义双流版的 result / chunk / aggregate 结构

---

## 4. 第一版推荐范围

第一版建议刻意收敛，不要一上来把所有扫描轴都打开。

推荐第一版只支持：

- 扫 `Eb/N0`
- `kChunkNumInfoBits`
- `kTargetPostErrors`
- `kMaxPostFecTotalBits`
- `kMaxTotalWorkers`
- `kMaxInflightChunksPerPoint`

其它参数全部固定，并且默认值直接与
[apps/ofec_two_stream_shared_decoder.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_decoder.cpp)
一致，包括：

- `kDecoderName`
- `kChaseL`
- `kChaseNTestOverride`
- `kChaseTopkKeep`
- `kChaseGroupMinimaBits`
- `kAlphaExplicit`
- `kBetaExplicit`
- `kSisoActiveList`
- `kHybridEnable` 及其相关配置
- `kMuxSchedulingMode`
- `kMuxPriorityRule`
- `kEnableEarlyStop`

这样做的好处是：

- one-shot app 和 sweep app 很容易对齐
- 如果扫出来的某个 `Eb/N0` 点结果异常，可以马上回到单点 app 复现
- 第一版问题更容易定位，不会把“扫描框架问题”、“BER 曲线问题”和“内部可观测性分析”混在一起

---

## 5. 顶层实现结构

`ofec_two_stream_shared_sweep.cpp` 建议沿用 `ofec_sweep3.cpp` 的文件组织方式：

1. 顶部 include
2. `namespace {}` 里的用户可调常量区
3. 双流版 `ChunkResult / AggregatedPointResult / PointState`
4. seed 派生函数
5. 停止条件与聚合逻辑
6. `run_chunk(...)`
7. 全局调度器
8. CSV 输出
9. `main()`

建议的骨架如下：

```cpp
#include <...>

#include "mux_bypass_edges.hpp"
#include "newcode/ofec_sweep_runner.hpp"
#include "newcode/two_stream_shared_runner.hpp"
#include "newcode/utils/now_stamp.hpp"
#include "ofec_sweep_detail.hpp"

namespace {

// ======== 用户可调参数区域 ========
// 与 ofec_two_stream_shared_decoder.cpp 对齐的固定参数
// 再加上 sweep3 风格的 chunk / stop / worker 参数

enum class StopReason { ... };

struct ChunkResult { ... };
struct AggregatedStreamStats { ... };
struct AggregatedPointResult { ... };
struct PointState { ... };
struct ActiveChunk { ... };

int derive_seed(...);
ChunkResult run_chunk(...);
void accumulate_chunk(...);
bool mark_stop_if_needed(...);
std::vector<AggregatedPointResult> run_low_ber_scenarios_global(...);

}  // namespace

int main() {
  ...
}
```

这里最重要的一点是：

> 主流程解码本身不在这个 app 里实现。  
> app 只负责“装配配置 + 调度 chunk + 聚合输出”。

---

## 6. 结果结构建议

## 6.1 `ChunkResult`

每个 chunk 的结果建议定义为：

```cpp
struct ChunkResult {
  std::size_t chunk_index = 0;
  int bitgen_seed_a = 0;
  int channel_seed_a = 0;
  int bitgen_seed_b = 0;
  int channel_seed_b = 0;
  newcode::two_stream_shared::Result result;
};
```

这里和单流 `sweep3` 最大的区别是：

- 需要记录四个 seed
- `result` 不再是 `newcode::PipelineResult`
- 而是 `newcode::two_stream_shared::Result`

## 6.2 `AggregatedStreamStats`

建议单独抽一个“双流公共统计结构”，避免 A/B 两边重复写三套字段：

```cpp
struct AggregatedStreamStats {
  newcode::BerStats pre_fec{};
  newcode::BerStats pre_fec_quantized_hard{};
  newcode::BerStats post_fec{};
  bool has_pre_fec_quantized_hard = false;
};
```

## 6.3 `AggregatedPointResult`

每个扫描点的聚合结果建议定义为：

```cpp
struct AggregatedPointResult {
  ofec_sweep::detail::SweepScenario scenario;
  AggregatedStreamStats stream_a;
  AggregatedStreamStats stream_b;
  std::size_t chunks_completed = 0;
  StopReason stop_reason_a = StopReason::NoData;
  StopReason stop_reason_b = StopReason::NoData;
  double elapsed_seconds = 0.0;
};
```

这里保留 `scenario`，是为了继续沿用现有 sweep 的场景命名和 CSV 记录方式。

---

## 7. 停止条件建议

## 7.1 推荐原则

双流时建议按“两个流都满足停止条件”才结束该点。

也就是说：

- A 满足了，不代表整个点结束
- B 满足了，也不代表整个点结束
- 只有 A 和 B 都满足停止条件，或者 hit 到全局上限，才停止

这样可以保证：

- A/B 两路的 BER 估计都足够稳定
- 不会出现其中一路统计明显偏薄的情况

## 7.2 推荐实现口径

建议先复用 `sweep3` 的三类停止依据：

1. `TargetPostErrors`
2. `MaxPostFecBits`
3. `ZeroErrorUpperBound`

但判定方式改为“对 A/B 分别判定，再取全点完成”：

- `stop_reason_a = evaluate_stop_reason(stream_a)`
- `stop_reason_b = evaluate_stop_reason(stream_b)`

全点停止条件建议是：

```text
done(point) = (stop_reason_a != NoData) && (stop_reason_b != NoData)
```

这里的一个细节是：

- A/B 两路每个 chunk 的 `post_fec.total` 通常相同
- 因此 `MaxPostFecBits` 很多时候会同时触发
- 但 `TargetPostErrors` 不一定同时触发

所以分路判断是必要的。

## 7.3 第一版不建议做的复杂逻辑

第一版先不要引入下面这些更复杂的策略：

- 按 A/B 中更差一路来动态加权 worker
- 对 A/B 单独分配不同 chunk 大小
- 某一路满足后只统计另一侧

这些后面都能加，但不是第一版要解决的问题。

---

## 8. 场景与 seed 设计

## 8.1 场景展开仍建议复用单流 `SweepScenario`

第一版建议仍然用：

- `ofec_sweep::SweepParameterConfig`
- `ofec_sweep::detail::build_scenarios(...)`

原因是 two-stream sweep 第一版主要只扫 `Eb/N0`，而这些公共工具已经够用了。

## 8.2 双流 seed 的推荐口径

需要额外补的是 two-stream 自己的四个基础 seed：

- `kBitgenSeedA`
- `kChannelSeedA`
- `kBitgenSeedB`
- `kChannelSeedB`

然后对每个 chunk 分别派生：

- `bitgen_seed_a(chunk)`
- `channel_seed_a(chunk)`
- `bitgen_seed_b(chunk)`
- `channel_seed_b(chunk)`

推荐继续复用 `sweep3` 的 `derive_seed(...)` 混合方式，只是对四路 seed 分别用不同 `salt`。

例如：

- `0x13579bdfU` 给 A bitgen
- `0x2468ace0U` 给 A channel
- `0x10293847U` 给 B bitgen
- `0x56473829U` 给 B channel

这样有两个好处：

- 每个 chunk 之间仍然独立
- A/B 两路不会因为偷懒共用 seed 而引入额外相关性

---

## 9. `run_chunk(...)` 的建议实现

`run_chunk(...)` 是整个 new app 的关键替换点。

单流 `sweep3` 当前做的是：

1. 从 `scenario` 构造 `Params`
2. 设置 `NUM_INFO_BITS`
3. 设置 bitgen/channel seed
4. 构造 `PipelineConfig`
5. 调 `newcode::run_pipeline(...)`

two-stream 版本建议改为：

1. 从固定常量区装配 `ofec_single::Config`
2. 通过现有 `build_params(...)` 和 `build_pipeline_config(...)` 得到基础 `Params + PipelineConfig`
3. 用 `scenario` 覆盖该点的 `Eb/N0` 和需要扫描的公共参数
4. 把 `NUM_INFO_BITS` 改为 `kChunkNumInfoBits`
5. 给 A/B 分别设置本 chunk 的 4 个 seed
6. 组装 `newcode::two_stream_shared::Config`
7. 调 `newcode::two_stream_shared::run_two_stream_shared(...)`

建议明确做到这一点：

> `ofec_two_stream_shared_sweep.cpp` 与 `ofec_two_stream_shared_decoder.cpp`
> 尽量使用同一套参数装配路径。

这样 sweep 跑出来的点，才能方便回切到单点 app 复现。

---

## 10. 聚合逻辑建议

每个 chunk 完成后，聚合器主要需要做两类事情。

## 10.1 聚合 BER 统计

分别对 A/B 做：

- `pre_fec.errors += ...`
- `pre_fec.total += ...`
- `post_fec.errors += ...`
- `post_fec.total += ...`

然后实时更新 `ber`。

## 10.2 更新停止状态

每次积完一个 chunk 后：

1. 分别检查 A 和 B 的停止原因
2. 若两边都已有停止原因，则停止继续 launch 新 chunk
3. 等 inflight chunk 全部收完后，完成该点

这个逻辑和 `sweep3` 的结构是一样的，只是从“单个 stop_reason”改成了“双 stop_reason”。

---

## 11. CSV 输出建议

CSV 的主目标应该是支撑 BER 曲线分析。

推荐至少包含这些字段：

- `timestamp`
- `run_id`
- `scenario`
- `decoder_name`
- `ebn0_db`

- `pre_ber_a`
- `pre_errs_a`
- `pre_total_a`
- `post_ber_a`
- `post_errs_a`
- `post_total_a`

- `pre_ber_b`
- `pre_errs_b`
- `pre_total_b`
- `post_ber_b`
- `post_errs_b`
- `post_total_b`

- `chunk_num_info_bits`
- `chunks_completed`
- `stop_reason_a`
- `stop_reason_b`
- `elapsed_seconds`

- 当前点的关键参数快照
  - `alpha_list`
  - `beta_list`
  - `siso_active_list`
  - `chase_L`
  - `chase_n_test`
  - `mux_scheduling_mode`
  - `mux_priority_rule`
  - `hybrid_enable`
  - `hybrid_classifier_mode`
  - `hybrid_siso_backfill_mode`

如果第一版不想一上来写太长的 CSV header，也至少要保证：

- BER
- stop reason
- 参数快照

这三大类一定在。

---

## 12. 第一版推荐落地步骤

## 12.1 第一步：先搭最小骨架

先新增：

- [apps/ofec_two_stream_shared_sweep.cpp](/home/zsr71/projects/newcode_two_stream_shared_chase/apps/ofec_two_stream_shared_sweep.cpp)

并完成这些最小能力：

- 固定参数区
- `Eb/N0` 扫描
- chunk 调度
- A/B BER 聚合
- CSV 输出

先不扩展太多扫描轴。

## 12.2 第二步：补更多扫描轴

第一批建议打开的扫描轴是：

- `decoder_name_candidates`
- `chase_l_candidates`
- `mux_scheduling_mode_candidates`
- `mux_early_stop_priority_rule_candidates`

不建议第一批就把所有 early-stop/hybrid 轴都一起开，否则 scenario 数会迅速膨胀。

---

## 13. 关于命名的说明

从长期维护角度看，`ofec_two_stream_shared_sweep.cpp` 比 `ofec_sweep4.cpp` 这种名字更清楚，因为它直接表明了：

- 这是 `two_stream`
- 这是 `shared`
- 这是 `sweep`

因此即使它现在仍然偏实验性质，这个命名也是合理的。

后续如果它逐渐成为主入口之一，这个名字也不需要再改。

---

## 14. 总结

`ofec_two_stream_shared_sweep.cpp` 的最佳定位不是“再造一个新的解码器”，而是：

- 用 `ofec_two_stream_shared_decoder.cpp` 提供双流共享主流程
- 用 `ofec_sweep3.cpp` 提供 low-BER chunk 扫描框架
- 在两者之间补一层双流版的参数装配和结果聚合

第一版最推荐的实现策略是：

- 只扫 `Eb/N0`
- 其它参数全部与单点 two-stream app 对齐
- 先把双流 BER 曲线稳定跑通

等这个版本稳定后，再逐步扩展到：

- 更多扫描轴
- 更复杂的自动停止与后处理逻辑

这样推进，工程风险最低，也最方便和现有单点 two-stream 验证入口相互对照。
