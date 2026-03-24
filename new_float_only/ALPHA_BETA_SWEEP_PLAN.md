# new_float_only `alpha/beta` 扫描实现计划

## 目标

在 [`apps`](/home/zsr71/projects/newcode/new_float_only/apps) 下新增一个独立的 `cpp` 入口，用于扫描 `ALPHA_LIST` / `beta_list` 参数，功能上参考老版本的 [`ofec_sweep2.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep2.cpp)，但实现上尽量复用 `new_float_only` 已有的 `run_single()` / `run_seed_sweep()` 能力，而不是再复制一套完整 pipeline。

预期新增目标类似：

- `apps/ofec_alpha_beta_sweep_float.cpp`
- 可执行文件名例如 `ofec_alpha_beta_sweep_float`

## 总体思路

旧版 `ofec_sweep2` 的核心职责可以拆成 4 部分：

1. 生成 `alpha/beta` 候选序列
2. 把每组候选转成一次可运行的 decoder 配置
3. 跑实验并记录 BER
4. 做 CSV/日志汇总，必要时分阶段筛选

在 `new_float_only` 里，最合理的迁移方式不是把旧文件直接搬过来，而是做一个“薄应用层”：

1. 在 `apps/` 中写扫描入口
2. 扫描入口内部构造一组 `SeedSweepConfig`
3. 每个候选调用现有的 [`run_seed_sweep()`](/home/zsr71/projects/newcode/new_float_only/include/new_float_only/seed_sweep_runner.hpp)
4. 把结果按 `pattern -> BER` 输出到新的 CSV

这样可以最大限度复用：

- 现有并行 seed sweep
- 现有 BER 聚合
- 现有日志/CSV 路径管理
- 现有 `Params` 校验逻辑

## 和旧版 `ofec_sweep2` 的差异

`new_float_only` 是 float-only、plain-only 的精简版本，因此实现上会主动做几处收缩：

1. 不保留量化相关参数
2. 不保留 decoder 名称切换
3. 不保留 interleaver 名称切换
4. 默认只扫 `float plain` 路径
5. 优先复用 `SeedSweepConfig`，不单独再造一套 sweep runner

也就是说，新版会保留“扫 `alpha/beta` 网格”的能力，但不会把老版的所有实验维度一起迁过来。

## 第一版建议功能

第一版先做一个够用的版本，不追求一次把旧版所有花样全部搬完。

建议第一版具备：

1. 固定 `CHASE_L`
2. 固定 `Eb/N0`
3. 固定 `trial_count`
4. 固定 `bits_per_symbol`
5. 通过常量定义一组 `alpha_low/high`、`beta_low/high` 网格
6. 自动生成每个 tile 的 `ALPHA_LIST` / `beta_list`
7. 对每组参数调用一次 `run_seed_sweep()`
8. 输出 summary CSV
9. 在控制台打印当前 pattern 和最终 BER

第一版先不做：

1. 两阶段筛选
2. keep-ratio 裁剪
3. 复杂 gamma 形状扫描
4. 参数组合并行调度
5. 命令行参数解析

原因很简单：先把“能稳定扫一批参数并写结果”做出来，再决定要不要上复杂调度。

## 推荐的数据结构

建议在新应用文件里定义两个轻量结构：

```cpp
struct Shape {
  float alpha_low;
  float alpha_high;
  float beta_low;
  float beta_high;
};

struct Evaluation {
  Shape shape;
  std::vector<float> alpha_list;
  std::vector<float> beta_list;
  new_float_only::SeedSweepResult result;
};
```

用途：

- `Shape` 表示一组待扫描的参数边界
- `Evaluation` 保存“输入参数 + 扫描结果”

如果后面要加入 `gamma_alpha/gamma_beta` 或多阶段筛选，再在 `Shape` 里扩展即可。

## `alpha/beta` 序列生成

旧版 `ofec_sweep2` 通过：

- `low`
- `high`
- `gamma`
- `tile_count`

来生成每个 tile 的参数序列。

`new_float_only` 第一版建议先保留这个建模方式，但默认把 `gamma=1.0` 写死，先只做线性插值：

```cpp
std::vector<float> build_sequence(float low, float high, size_t count);
```

例如 `TILES_PER_WIN = 4` 时：

- `alpha_low = 0.2`
- `alpha_high = 0.8`

生成：

```cpp
{0.2, 0.4, 0.6, 0.8}
```

这样和现在 [`ofec_seed_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_seed_sweep_float.cpp) 里写死的常量风格一致，迁移成本最低。

## 配置复用方式

每次评估一个候选时，流程大致如下：

1. 构造基础 `SeedSweepConfig`
2. 填入公共配置
3. 覆盖当前候选的 `decoder.ALPHA_LIST`
4. 覆盖当前候选的 `decoder.beta_list`
5. 调用 `run_seed_sweep(config)`

伪代码如下：

```cpp
new_float_only::SeedSweepConfig config = build_base_config();
config.decoder.ALPHA_LIST = alpha_list;
config.decoder.beta_list = beta_list;
const auto result = new_float_only::run_seed_sweep(config);
```

这里 `build_base_config()` 会集中设置：

- `label`
- `ebn0_db`
- `trial_count`
- `bitgen_seed_base`
- `channel_seed_base`
- `normalize_extrinsic`
- `CHASE_L`
- `CHASE_NTEST`
- `HARD_TILE_LIST`
- `debug_trace`

这样应用层只改扫描参数，不改基础实验框架。

## 输出设计

建议第一版生成一个 CSV，字段至少包括：

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

- `alpha_list` / `beta_list` 建议直接写成带分号的字符串
- 文件路径建议放在 `data/` 下
- 文件名建议带时间戳，避免覆盖

例如：

```text
data/ofec_alpha_beta_sweep_float_<timestamp>.csv
```

## 并发策略

第一版建议：

- 参数模式之间串行
- 单个模式内部复用 `run_seed_sweep()` 自己的并行 trial 机制

原因：

1. 现有 `run_seed_sweep()` 已经有 worker 限流
2. 如果“pattern 并行”再叠加“trial 并行”，很容易线程打爆
3. 串行 pattern 更容易看日志和定位异常

如果后面确认 CPU 还有富余，再考虑加一层 pattern 级并行。

## 文件改动计划

预计会改这些地方：

1. 新增 [`apps/ofec_alpha_beta_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_alpha_beta_sweep_float.cpp)
2. 修改 [`CMakeLists.txt`](/home/zsr71/projects/newcode/new_float_only/CMakeLists.txt)

第一版尽量不改：

1. `seed_sweep_runner.cpp`
2. `pipeline_runner.cpp`
3. `ofec_decoder` 主流程

只有在发现应用层很难复用现有接口时，才考虑把一小部分公共逻辑抽到 `src/seed_sweep/`。

## 可能的第二阶段增强

如果第一版跑通，后面再加这些能力：

1. 两阶段筛选
   第一阶段少量 bit 快扫
   第二阶段只保留前若干模式精扫

2. `gamma_alpha / gamma_beta`
   允许非线性地生成 tile 序列

3. 显式 pattern label
   让 CSV 和日志里更容易定位参数组合

4. pattern 级并行
   用受控线程池并发多组参数

5. 命令行参数
   不再只能改源码常量

## 我准备按这个顺序实现

1. 先新增一个最小可运行版 `apps/ofec_alpha_beta_sweep_float.cpp`
2. 只支持线性 `alpha/beta` 序列
3. 只输出一个 summary CSV
4. 接到 `CMakeLists.txt`
5. 编译并跑通一轮
6. 再看是否需要补两阶段筛选或更复杂 pattern 生成

## 当前判断

最稳妥的方案是：

- 不复制旧版 `ofec_sweep2` 的整套 runner
- 只迁移“pattern 生成 + 结果汇总”这层应用逻辑
- 底层实验执行统一复用 `new_float_only::run_seed_sweep()`

这样代码量最小，和 `new_float_only` 当前结构也最一致。
