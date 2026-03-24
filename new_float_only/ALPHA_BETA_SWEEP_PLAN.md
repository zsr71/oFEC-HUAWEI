# `new_float_only` `alpha/beta` 扫描实现计划

## 目标

在 [`apps`](/home/zsr71/projects/newcode/new_float_only/apps) 下新增一个独立入口，用于扫描 `ALPHA_LIST` / `beta_list` 参数，功能上参考老版本的 [`ofec_sweep2.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep2.cpp)，但实现范围严格限制在 [`new_float_only`](/home/zsr71/projects/newcode/new_float_only) 内。

预期新增目标类似：

- `apps/ofec_alpha_beta_sweep_float.cpp`
- 可执行文件名例如 `ofec_alpha_beta_sweep_float`

## 总体思路

旧版 `ofec_sweep2` 的核心职责可以拆成 4 部分：

1. 生成 `alpha/beta` 候选序列
2. 把每组候选转成一次可运行的 decoder 配置
3. 跑实验并记录 BER
4. 做 CSV/日志汇总，必要时分阶段筛选

在 `new_float_only` 里，不建议把旧文件直接搬过来，而是做一个更小的统一框架：

1. 在 `apps/` 中写扫描入口
2. 先生成一份固定共享的 seed schedule
3. 再生成一组参数 pattern
4. 把 `pattern × seed` 展开成 task 列表
5. 并行执行这些 task
6. 最后按 `pattern` 聚合 BER，并输出 CSV

这样既能保证不同参数模式使用相同 seed，又能把 `seed sweep` 和 `parameter sweep` 收敛到同一套底层调度模型。

## 和旧版 `ofec_sweep2` 的差异

`new_float_only` 是 float-only、plain-only 的精简版本，因此实现上主动收缩范围：

1. 不保留量化相关参数
2. 不保留 decoder 名称切换
3. 不保留 interleaver 名称切换
4. 默认只扫 `float plain` 路径
5. 扫描层自己控制参数模式、seed 调度和 task 并行

也就是说，新版只保留“扫 `alpha/beta` 网格”的能力，不把旧版所有实验维度一起迁进来。

## 统一到 Task Scheduler

后续想统一“seed sweep”和“parameter sweep”，底层不应该继续区分“这是一批 seed”还是“这是一批参数”，而应该统一成“很多个单次实验 task”。

建议抽象为：

```cpp
struct SweepPattern {
  std::size_t pattern_index;
  std::string pattern_label;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
};

struct SweepSeedPair {
  std::size_t seed_index;
  int bitgen_seed;
  int channel_seed;
};

struct SweepTask {
  std::size_t task_index;
  std::size_t pattern_index;
  std::size_t seed_index;
  SweepPattern pattern;
  SweepSeedPair seed_pair;
};
```

这个模型下：

- 纯 `seed sweep` = `1 pattern × N seeds`
- 纯 `parameter sweep` = `N patterns × M seeds`
- 将来要加 `Eb/N0 sweep` 也可以继续扩展 task 维度

## 为什么最小执行单元要用 `run_pipeline()`

文档里不再建议并发直接调用 `run_single()`，原因很明确：

1. [`run_single()`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/single_runner.hpp) 会创建日志文件
2. 当前日志文件名依赖秒级时间戳
3. 多个并发 task 在同一秒启动时，存在文件名冲突风险

所以统一 `task scheduler` 时，底层更适合直接复用：

- [`run_pipeline()`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/pipeline_runner.hpp)

由 task 执行器自行构造：

- `Params`
- `PipelineConfig`
- `label`
- `ebn0_db`

这样既能避开并发日志冲突，也能让 `seed sweep` 和 `parameter sweep` 共用完全一致的最小执行单元。

## 第一版建议功能

第一版先做一个够用的版本，不追求一次搬完旧版所有能力。

建议具备：

1. 固定 `CHASE_L`
2. 固定 `Eb/N0`
3. 固定 `trial_count`
4. 固定 `bits_per_symbol`
5. 通过常量定义一组 `alpha_low/high`、`beta_low/high` 网格
6. 自动生成每个 tile 的 `ALPHA_LIST` / `beta_list`
7. 预先生成一份固定共享的 seed 列表
8. 把所有 `(pattern, seed)` 组合展开成 task 列表
9. 并行执行这些 task
10. 输出 summary CSV
11. 在控制台打印当前 pattern 和最终 BER

第一版先不做：

1. 两阶段筛选
2. keep-ratio 裁剪
3. 复杂 gamma 形状扫描
4. 命令行参数解析

原因很简单：先把“能稳定扫一批参数并写结果”做出来，再决定要不要上复杂调度。

## 推荐的数据结构

扫描 app 可以先定义轻量 pattern 结构：

```cpp
struct Shape {
  float alpha_low;
  float alpha_high;
  float beta_low;
  float beta_high;
};
```

公共 task scheduler 则负责：

```cpp
struct SweepTaskResult {
  std::size_t task_index;
  std::size_t pattern_index;
  std::size_t seed_index;
  int bitgen_seed;
  int channel_seed;
  new_float_only::BerStats pre_fec;
  new_float_only::BerStats post_fec;
};
```

按职责分层：

- app 层负责定义 pattern 空间
- task scheduler 负责执行 `(pattern, seed)` task
- 聚合层负责把 task result 汇总成 pattern 级 BER

## `alpha/beta` 序列生成

旧版 `ofec_sweep2` 通过：

- `low`
- `high`
- `gamma`
- `tile_count`

来生成每个 tile 的参数序列。

`new_float_only` 第一版建议先保留这个建模方式，但默认把 `gamma=1.0` 写死，先只做线性插值：

```cpp
std::vector<float> build_sequence(float low, float high, std::size_t count);
```

例如 `TILES_PER_WIN = 4` 时：

- `alpha_low = 0.2`
- `alpha_high = 0.8`

生成：

```cpp
{0.2, 0.4, 0.6, 0.8}
```

这样和现有 [`ofec_seed_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_seed_sweep_float.cpp) 里写死常量的风格一致，迁移成本最低。

## 配置复用方式

每次评估一个候选时，流程大致如下：

1. 构造基础 decoder 配置
2. 补齐公共链路配置
3. 用 pattern 覆盖 `decoder.ALPHA_LIST`
4. 用 pattern 覆盖 `decoder.beta_list`
5. 用固定 seed schedule 生成 task
6. 并行执行 task
7. 按 `pattern_index` 聚合 BER

伪代码如下：

```cpp
auto seeds = build_shared_seed_schedule();
auto tasks = build_tasks(patterns, seeds);
auto task_results = run_tasks_parallel(tasks);
auto summary = aggregate_by_pattern(task_results);
```

进入 task 内部后，再构造：

```cpp
Params params = base_decoder;
params.ALPHA_LIST = task.pattern.alpha_list;
params.beta_list = task.pattern.beta_list;

PipelineConfig pipeline;
pipeline.bitgen_seed = task.seed_pair.bitgen_seed;
pipeline.channel_seed = task.seed_pair.channel_seed;

auto result = run_pipeline(params, pipeline, task_label, ebn0_db);
```

## 输出设计

建议第一版生成一个 summary CSV，字段至少包括：

1. `run_id`
2. `pattern_label`
3. `ebn0_db`
4. `trial_count`
5. `alpha_list`
6. `beta_list`
7. `post_ber`
8. `post_errors`
9. `post_total`
10. `pre_ber`
11. `pre_errors`
12. `pre_total`

其中：

- `alpha_list` / `beta_list` 直接写成带分号的字符串
- 文件路径放在 `data/`
- 文件名带时间戳，避免覆盖

例如：

```text
data/ofec_alpha_beta_sweep_float_<timestamp>.csv
```

## 并发策略

并发单元定义为一个 `(pattern, seed)` task。

执行方式：

1. 先生成一份固定共享的 seed schedule
2. 再生成参数 pattern 列表
3. 展开成 `pattern × seed` 的 task 列表
4. 线程池直接并行执行这些 task
5. 所有 task 完成后，再按 `pattern` 聚合结果

例子：

- 固定共享 seed 为 3 对：
  - `(bitgen=1001, channel=2001)`
  - `(bitgen=1002, channel=2002)`
  - `(bitgen=1003, channel=2003)`
- 要比较 4 组参数模式：A/B/C/D

正确执行方式是：

1. A、B、C、D 都使用同样这 3 对 seed
2. 展开后总共有 `4 × 3 = 12` 个 task
3. 线程池直接并行执行这 12 个 task
4. 最后把 `(A,S1)(A,S2)(A,S3)` 聚合成 A 的 BER，B/C/D 同理

这样 A/B/C/D 的 BER 差异主要来自参数本身，而不是来自不同 seed 采样。

## 文件改动计划

预计会改这些地方：

1. 新增 [`include/new_float_only/sweep_task_runner.hpp`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/sweep_task_runner.hpp)
2. 新增 `src/seed_sweep/sweep_task_runner.cpp`
3. 修改 [`src/seed_sweep/seed_sweep_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/seed_sweep/seed_sweep_runner.cpp)
4. 新增 [`apps/ofec_alpha_beta_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_alpha_beta_sweep_float.cpp)
5. 修改 [`CMakeLists.txt`](/home/zsr71/projects/newcode/new_float_only/CMakeLists.txt)

改动范围仍然只限于 [`new_float_only`](/home/zsr71/projects/newcode/new_float_only)。

## 可能的第二阶段增强

如果第一版跑通，后面再加这些能力：

1. 两阶段筛选
   第一阶段少量 bit 快扫，第二阶段只保留前若干模式精扫
2. `gamma_alpha / gamma_beta`
   允许非线性地生成 tile 序列
3. 显式 pattern label 模板
   让 CSV 和控制台输出更容易定位参数组合
4. 多维联合扫描
   例如同时扫描 `Eb/N0` 和 `alpha/beta`
5. 命令行参数
   不再只能改源码常量

## 我准备按这个顺序实现

1. 先补统一 `task scheduler` 的公共头文件和实现
2. 让 `run_seed_sweep()` 先改为复用这套公共调度层
3. 再新增最小可运行版 `apps/ofec_alpha_beta_sweep_float.cpp`
4. 只支持线性 `alpha/beta` 序列
5. 输出一个 summary CSV
6. 接到 `CMakeLists.txt`
7. 编译并跑通一轮

## 当前判断

最稳妥的方案是：

- 不复制旧版 `ofec_sweep2` 的整套 runner
- 先把底层统一成 `task scheduler + run_pipeline()`
- 再把 `seed sweep` 和 `parameter sweep` 都建立在这个公共底座上

这样代码量可控，后面继续扩展也不会出现两套 sweep 逻辑各自演化的问题。
