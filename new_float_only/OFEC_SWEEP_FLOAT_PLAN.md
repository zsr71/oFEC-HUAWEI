# `new_float_only` 简化版 Eb/N0 Sweep App 实现方案

## 目标

在 [`new_float_only/apps`](/home/zsr71/projects/newcode/new_float_only/apps) 下新增一个独立 app，功能收敛为：

- 给定一组固定的 float-only decoder / 链路参数
- 扫描一段 `Eb/N0`
- 并行跑任务
- 汇总每个 `Eb/N0` 点的 BER 结果

它的风格参考主仓 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)，但能力范围明显更小，只使用 [`new_float_only`](/home/zsr71/projects/newcode/new_float_only) 现有代码，不依赖主仓 `newcode` 的 sweep runner。

建议新增：

- `new_float_only/apps/ofec_ebn0_sweep_float.cpp`

可执行文件名建议：

- `ofec_ebn0_sweep_float`

## 功能边界

这个新 app 只做一件事：

- 固定参数
- 展开多个 `Eb/N0` 点
- 对每个 `Eb/N0` 点跑一批 seed trial
- 并行执行
- 输出每个 `Eb/N0` 点的聚合 BER

明确不做：

- 不做 `alpha/beta` 网格搜索
- 不做 `CHASE_L` / `CHASE_NTEST` 候选笛卡尔积
- 不做 decoder 变体切换
- 不做 early-stop / MUX 多轴扫描
- 不做“best scenario”搜索

也就是说，这个 app 不是“统一多轴 sweep”，而是“固定参数下的 `Eb/N0` sweep”。

## 为什么这样收敛

你现在要的主功能其实很明确：

1. 先定一组参数
2. 看这组参数在不同 `Eb/N0` 下的表现
3. 为了效率，把这些任务并行跑掉

在 `new_float_only` 里，这个需求不需要做一个像主仓 `ofec_sweep.cpp` 那样的大型 scenario 系统。  
更合理的做法是做一个轻量 runner：

- `Eb/N0` 是唯一 sweep 轴
- seed 只是 Monte Carlo 重复试验维度
- decoder 参数全部固定

这样更简单，也更贴合 `new_float_only` 当前已有基础设施。

## 可复用的现有基础设施

### 1. float-only 参数和单次链路

- [`new_float_only/include/new_float_only/params.hpp`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/params.hpp)
- [`new_float_only/include/new_float_only/pipeline_runner.hpp`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/pipeline_runner.hpp)

这里已经能承载固定参数：

- `DecoderConfig`
- `Params`
- `run_pipeline()`

### 2. 通用并行 task scheduler

- [`new_float_only/include/new_float_only/sweep_task_runner.hpp`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/sweep_task_runner.hpp)
- [`new_float_only/src/seed_sweep/sweep_task_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/seed_sweep/sweep_task_runner.cpp)

它已经能做：

- 构造 seed schedule
- 把任务列表并行跑掉
- 聚合 BER

虽然它现在的抽象名叫 `SweepPattern`，但完全可以把“每个 `Eb/N0` 点”当成一个 pattern 来复用。

## 推荐实现思路

最务实的实现不是重新造线程调度，而是包一层新的 runner，把 `Eb/N0` 这一轴映射到现有 task scheduler。

建议新增：

- [`new_float_only/include/new_float_only/ebn0_sweep_runner.hpp`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/ebn0_sweep_runner.hpp)
- `new_float_only/src/seed_sweep/ebn0_sweep_runner.cpp`

app 入口：

- [`new_float_only/apps/ofec_ebn0_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_ebn0_sweep_float.cpp)

## 推荐的数据模型

### 配置结构

```cpp
struct Ebn0SweepConfig {
  std::string label = "ofec_ebn0_sweep_float";

  float ebn0_start = 3.0f;
  float ebn0_end = 3.0f;
  int ebn0_points = 1;

  unsigned bits_per_symbol = 1;
  bool generate_random_bits = true;
  bool normalize_extrinsic = true;
  DecoderConfig decoder{};

  std::size_t trial_count = 1;
  int bitgen_seed_base = 56456;
  int channel_seed_base = 57112;
  std::vector<int> bitgen_seeds;
  std::vector<int> channel_seeds;

  unsigned max_workers_override = 0;
  bool quiet_pipeline = true;
  bool quiet_logs = false;
  bool write_summary_csv = true;
  bool write_trial_csv = true;
};
```

### 单点汇总结果

```cpp
struct Ebn0SweepPointResult {
  std::size_t ebn0_index = 0;
  float ebn0_db = 0.0f;
  std::size_t trials_completed = 0;
  BerStats pre_fec;
  BerStats post_fec;
  std::size_t pre_frame_error_trials = 0;
  std::size_t post_frame_error_trials = 0;
};
```

### 总结果

```cpp
struct Ebn0SweepResult {
  std::vector<Ebn0SweepPointResult> points;
  std::string summary_csv_path;
  std::string trial_csv_path;
};
```

## 任务展开方式

建议把每个 `Eb/N0` 点看成一个 pattern，再和 seed schedule 做笛卡尔积：

- `Eb/N0 points × seed trials`

例如：

- `Eb/N0 = {3.0, 3.1, 3.2}`
- `trial_count = 4`

则总任务数是：

- `3 × 4 = 12`

每个 task 对应：

- 一个固定 `Eb/N0`
- 一对固定 `(bitgen_seed, channel_seed)`
- 一组固定 decoder 参数

这正好符合现有 `run_sweep_tasks()` 的并行模型。

## 复用现有 scheduler 的办法

`new_float_only` 现有 scheduler 里，`pattern` 主要携带：

- `pattern_index`
- `pattern_label`
- `alpha_list`
- `beta_list`

对这个新 app 来说，decoder 参数是固定的，所以：

- `alpha_list` / `beta_list` 直接沿用固定 decoder 里的值
- `pattern_label` 可以写成 `EbN0_3.10`
- `pattern_index` 对应 `ebn0_index`

真正的 `Eb/N0` 数值则在 runner 里单独维护。  
如果你觉得复用现有 `SweepPattern` 语义太绕，也可以新建一套 `Ebn0Task` 结构，但第一版没必要。

## 推荐 runner 流程

建议 `run_ebn0_sweep(config)` 的流程如下：

1. 规范化 `DecoderConfig`
2. 构造 `Eb/N0` 网格
3. 构造共享 seed schedule
4. 展开成 `(ebn0, seed)` 任务
5. 并行执行任务
6. 按 `Eb/N0` 聚合结果
7. 输出 summary/trial CSV
8. 返回结果

伪代码如下：

```cpp
Params params = normalize_sweep_decoder_config(config.decoder);
auto ebn0_values = build_ebn0_grid(config.ebn0_start, config.ebn0_end, config.ebn0_points);
auto seed_schedule = build_sweep_seed_schedule(...);
auto tasks = build_ebn0_tasks(ebn0_values, seed_schedule);
auto task_results = run_ebn0_tasks(config, params, tasks);
auto summary = aggregate_by_ebn0(task_results);
write_csv(...);
return summary;
```

## app 入口建议结构

[`new_float_only/apps/ofec_ebn0_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_ebn0_sweep_float.cpp) 建议保持和现有 app 一样的风格：

```cpp
#include <iostream>

#include "new_float_only/ebn0_sweep_runner.hpp"

namespace {

constexpr const char* kLabel = "ofec_ebn0_sweep_float";
constexpr float kEbN0Start = 3.0f;
constexpr float kEbN0End = 3.4f;
constexpr int kEbN0Points = 9;
constexpr std::size_t kTrialCount = 4;

constexpr int kBitgenSeedBase = 1521867291;
constexpr int kChannelSeedBase = 998258255;
constexpr unsigned kBitsPerSymbol = 1;
constexpr bool kGenerateRandomBits = true;
constexpr bool kNormalizeExtrinsic = true;
constexpr unsigned kMaxWorkersOverride = 0;
constexpr bool kQuietPipeline = true;
constexpr bool kQuietLogs = false;

constexpr bool kNormalizeKnownPrefixTail = false;
constexpr int kChaseL = 6;
const std::vector<float> kAlphaList = {0.2f, 0.4f, 0.6f, 0.8f};
const std::vector<float> kBetaList = {0.2f, 0.4f, 0.6f, 0.8f};

}  // namespace

int main() {
  new_float_only::Ebn0SweepConfig config;
  // 填 config
  const auto result = new_float_only::run_ebn0_sweep(config);
  return 0;
}
```

重点是：

- 顶部常量区直接改参数
- `main()` 只做 config 装配
- 真正执行逻辑放到 runner 里

## 输出建议

建议输出两个 CSV。

### 1. summary CSV

每行一个 `Eb/N0` 点，字段至少包括：

- `run_id`
- `label`
- `ebn0_index`
- `ebn0_db`
- `trials_requested`
- `trials_completed`
- `pre_ber`
- `pre_errors`
- `pre_total`
- `pre_frame_error_trials`
- `post_ber`
- `post_errors`
- `post_total`
- `post_frame_error_trials`

### 2. trial CSV

每行一个 task，字段至少包括：

- `run_id`
- `ebn0_index`
- `ebn0_db`
- `trial_index`
- `bitgen_seed`
- `channel_seed`
- `pre_ber`
- `post_ber`

这样既能画 BER 曲线，也能排查单个 seed 异常点。

## 控制台输出建议

控制台不需要像主仓 `ofec_sweep.cpp` 那样打印“best scenario”，因为这里没有 scenario 竞争。  
更合适的输出是：

- 启动时打印 `Eb/N0` 范围、点数、trial 数、worker 数
- 运行中打印进度
- 完成后逐行打印每个 `Eb/N0` 的 post-FEC BER
- 最后打印 CSV 路径

## 不应支持的并发副作用

沿用现有 `sweep_task_runner` 的限制：

- 不支持 `DUMP_WORK_LLR`
- 不支持 `debug_trace.enable`

因为多个并发 task 会覆盖同一路径或混淆调试输出。

所以 runner 在进入并行执行前，应继续复用现有的校验逻辑。

## 最小实现步骤

建议按下面顺序落地：

1. 新增 `ebn0_sweep_runner.hpp`
2. 新增 `ebn0_sweep_runner.cpp`
3. 先实现 `build_ebn0_values()`
4. 先实现 `(ebn0, seed)` task 展开
5. 先实现并行执行和按 `Eb/N0` 聚合
6. 再补 summary/trial CSV
7. 最后新增 `apps/ofec_ebn0_sweep_float.cpp`

## 最小验收标准

做到下面这些，就可以认为第一版完成：

1. 能编译出 `ofec_ebn0_sweep_float`
2. 给定固定 decoder 参数后，能扫描一段 `Eb/N0`
3. 任务以 `(ebn0, seed)` 为单位并行执行
4. 能输出每个 `Eb/N0` 点的 pre/post BER
5. 能输出 summary CSV 和 trial CSV
6. 整个实现只使用 `new_float_only`

## 结论

这个新 app 不需要做成主仓 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 那种多轴大扫面器。  
对你现在的目标，最合适的方案是：

- 固定参数
- 只扫描 `Eb/N0`
- 并行跑 `(ebn0, seed)` 任务
- 输出 BER 曲线数据

这比做一个“大而全的 float-only sweep 框架”更直接，也更稳。
