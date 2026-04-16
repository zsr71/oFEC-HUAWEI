# 方案C为主的 BER 起始 Window 测试方案

## 1. 目标

当前 BER 统计会裁掉前若干个 window 作为 warm-up 区，并忽略最后一个 window。问题不在于“能不能裁”，而在于：

- 前面到底应该裁掉几个 window
- 从第几个 window 开始统计，结果才足够稳定
- 如果裁得过多，是否会浪费有效样本

本方案的目标是为“BER 从第几个 window 开始统计”提供实验依据，而不是先修改正式 BER 逻辑。

## 2. 总体思路

推荐以 **方案 C** 为主体，辅以 **方案 B** 和 **方案 A**：

- 方案 C：多 seed 重复后的 per-window 聚合统计
- 方案 B：从第 `k` 个 window 开始累计 BER 的稳定性分析
- 方案 A：逐 window BER 剖面观察

三者分工如下：

- 方案 A 负责看时域形态
- 方案 B 负责给出“从哪里开始算更合适”的判据
- 方案 C 负责把统计做稳，避免单次运行噪声过大

## 3. 为什么以方案 C 为主体

只看单次运行的逐 window BER，风险很明显：

- 单个 window 样本量有限
- 低 BER 区间下单 window 很容易出现 0 error
- 单 seed 结果容易受偶然性影响

方案 C 的核心优势是：对同一配置、同一 `window_idx`，跨多个 seed 聚合错误数和总比特数，再统一计算 BER：

```text
BER_window[w] = sum(errs_s[w]) / sum(bits_s[w])
```

这样得到的是统一口径下的“总错误数 / 总比特数”，比简单平均多个 BER 更合理，也更接近现有 runner 的 BER 聚合方式。

## 4. 测试对象与基本约束

第一版建议固定一组通信链路配置，只研究“BER 起始 window 选择”本身，不同时扫描太多变量。

建议固定的内容包括：

- 调制方式
- 译码器结构与参数
- early-stop 和 MUX 配置
- window 划分方式
- 现有 warm-up 机制保持不变

然后在若干代表性 `Eb/N0` 点上做测试，例如：

- 偏低 SNR 点：BER 较高
- 目标工作点：你真正关心的设计点
- 偏高 SNR 点：更接近稳态低 BER 区间

## 5. 推荐记录的数据

对每个 seed、每个 window，建议记录以下原始统计量：

- `window_idx`
- `pre_errs`
- `pre_bits`
- `post_errs`
- `post_bits`

第一版不必加过多辅助指标。只要 window 级错误数和总比特数能导出来，就足够支撑方案 A/B/C。

## 6. 三层分析结构

### 6.1 第一层：单 seed、逐 window 原始统计

这一层保留原始信息，用于做 trace 和异常定位。

对单个 seed，可定义：

```text
BER_s[w] = errs_s[w] / bits_s[w]
```

这层主要是中间数据，不作为最终结论。

### 6.2 第二层：多 seed 聚合的逐 window BER

对同一个 `window_idx`，在多个 seed 上聚合：

```text
BER_agg[w] = sum(errs_s[w]) / sum(bits_s[w])
```

这就是方案 C 的主体输出。它回答的问题是：

- 第 `w` 个 window 在统计意义上的误码表现如何

这一步最适合画“时域 BER 剖面图”。

### 6.3 第三层：从第 k 个 window 开始累计的 BER

对每个候选起点 `k`，定义：

```text
BER_ge_k = sum_{w>=k}(errs[w]) / sum_{w>=k}(bits[w])
```

这个量最关键，因为它直接对应工程问题：

- 如果从第 `k` 个 window 开始纳入 BER 统计，最终 BER 会是多少

## 7. 方案 A 的作用

方案 A 不是主要判据，但非常有价值。

逐 window BER 剖面可以帮助你观察：

- 前段是否存在明显的过渡区
- 某几个 window 是否异常高
- BER 是否单调收敛，还是存在回摆或局部波动

所以建议保留它作为辅助可视化，但不要只凭肉眼判断最终起点。

## 8. 推荐测试流程

### 步骤 1：固定配置

先只固定一组完整配置，避免多个变量同时变化。

### 步骤 2：选择若干代表性工作点

建议至少测 2 到 3 个 `Eb/N0` 点，避免结论只适用于单一点。

### 步骤 3：多 seed 重复运行

推荐分两轮：

- 第一轮：5 个 seed，快速看趋势
- 第二轮：10 到 20 个 seed，做正式判断

### 步骤 4：生成三类结果

输出三类结果：

- 逐 window 聚合 BER：`BER_agg[w]`
- 从第 `k` 个 window 开始累计 BER：`BER_ge_k`
- 相对参考 BER 的偏差：`delta(k)`

## 9. 推荐判据

### 9.1 参考 BER 的定义

需要先定义一个“稳态参考 BER”。第一版建议用下面两种方式之一：

- 后半段 windows 的累计 BER
- 从一个较保守起点开始的累计 BER，例如最后 50% windows

记为：

```text
BER_ref
```

### 9.2 相对偏差

对每个候选起点 `k`，定义：

```text
delta(k) = abs(BER_ge_k - BER_ref) / BER_ref
```

### 9.3 推荐规则

选择最小的 `k`，满足：

- `delta(k) < epsilon`，例如 5% 或 10%
- 从该 `k` 开始，后续若干个候选点也持续满足这一条件
- 从 `k` 开始累计的总比特数仍足够多

这个规则比“看图觉得平稳了”更适合后续固化到 app 中。

## 10. 建议输出图表

建议至少画三张图：

### 图 1：逐 window 聚合 BER

- 横轴：`window_idx`
- 纵轴：`BER_agg[w]`

作用：

- 观察时域上的 warm-up 区和稳态区

### 图 2：从第 k 个 window 开始累计 BER

- 横轴：`k`
- 纵轴：`BER_ge_k`

作用：

- 直接看裁掉前多少个 window 后，最终 BER 如何变化

### 图 3：相对参考 BER 偏差曲线

- 横轴：`k`
- 纵轴：`delta(k)`

作用：

- 直接找推荐起点

可选再补一张：

- 多 seed 单独曲线叠加图

## 11. 建议输出文件

建议分两类 CSV。

### 11.1 原始级别输出

每个 seed、每个 window 一条记录：

- `config_id`
- `ebn0_db`
- `seed`
- `window_idx`
- `pre_errs`
- `pre_bits`
- `post_errs`
- `post_bits`

### 11.2 聚合级别输出

按 window 输出：

- `window_idx`
- `agg_pre_errs`
- `agg_pre_bits`
- `agg_pre_ber`
- `agg_post_errs`
- `agg_post_bits`
- `agg_post_ber`

按候选起点 `k` 输出：

- `start_window_k`
- `post_errs_from_k`
- `post_bits_from_k`
- `post_ber_from_k`
- `reference_ber`
- `relative_deviation`
- `stable_flag`

最终建议输出：

- `recommended_start_window`
- `recommended_reason`
- `epsilon`
- `reference_definition`

## 12. 工程实现建议

### 12.1 第一版先做独立分析脚本

建议新增一个独立 app，而不是直接改现有 BER 逻辑。理由很简单：

- 不影响现有 runner
- 便于反复调整判据
- 输出更灵活

### 12.2 保持当前 warm-up 机制不变

第一阶段的目标是“先观察，再定规则”，不是立刻替换正式 warm-up 设置。

所以建议：

- 保持现有 warm-up 机制
- 只额外导出 per-window 统计
- 通过离线分析或半在线分析决定未来默认起点

### 12.3 脚本支持两种模式

建议脚本支持：

- 单 seed 调试模式
- 多 seed 统计模式

前者用于快速看现象，后者用于产出正式结论。

## 13. 可能遇到的问题

### 13.1 低 BER 区间下单 window 大量 0 error

这是正常现象。应对方式：

- 以多 seed 聚合为主
- 更重视 `BER_ge_k`
- 不要只看单个 `BER[w]`

### 13.2 后段样本量不足

如果总 window 数不多，裁掉前几个后，后面的统计量可能不够。

因此自动判据必须同时约束：

- `delta(k)`
- `bits_from_k`

### 13.3 不同工作点推荐起点不同

这是可能发生的。

处理方式有两种：

- 分工作点分别给建议
- 取一个更保守的统一值作为系统默认

工程上通常第二种更容易维护。

## 14. 推荐结论

针对当前问题，推荐的主方案是：

- 以 **方案 C** 为主体：多 seed 聚合逐 window BER
- 以 **方案 B** 为核心判据：从第 `k` 个 window 开始累计 BER
- 以 **方案 A** 为辅助观察：逐 window BER 剖面图

这样做的优势是：

- 既能看到现象
- 又能得到定量结论
- 还能顺着现有统计框架往下实现

## 15. 第一版落地建议

第一版先实现这些最关键的功能：

### 必做

- 多 seed、per-window 原始错误数统计
- 聚合后的 per-window BER
- `BER_ge_k` 曲线
- 自动推荐 `recommended_start_window`

### 选做

- 多 seed 单独曲线叠加
- 同时输出 pre-FEC / post-FEC
- 加入 frame error 作为辅助指标

### 暂不建议第一版做

- 复杂的变点检测算法
- 过于复杂的自适应判据
- 与主 runner 深度耦合的自动 warm-up 改写逻辑

## 16. 代码改造范围建议

这一部分的目标是明确：

- 哪些代码建议新增
- 哪些代码建议修改
- 每一块代码各自承担什么职责

总体原则是：

- **新增一个独立 app 做测试**
- **尽量复用现有 pipeline / BER 统计 / 配置接线**
- **不要直接改正式 BER 逻辑**

### 16.1 建议新增的代码

#### 1. 新增 app 入口

建议新增一个独立脚本，例如：

- `apps/ofec_ber_window_probe.cpp`

职责：

- 固定一组或若干组测试配置
- 控制多 seed 重复运行
- 调用现有 pipeline
- 收集每个 window 的原始统计
- 输出原始 CSV 和聚合 CSV
- 生成推荐的 `recommended_start_window`

这是第一版最核心的新增文件。

#### 2. 新增 window 级统计结构

建议新增一个轻量级的数据结构，例如放在：

- `include/newcode/pipeline_runner.hpp`
- 或新建一个专门的头文件，例如：
  - `include/newcode/ber_window_probe.hpp`

建议结构至少包含：

- `window_idx`
- `pre_errs`
- `pre_bits`
- `post_errs`
- `post_bits`

如果要支持多 seed 聚合，再补：

- `seed`
- `ebn0_db`
- `config_id`

#### 3. 新增分析与导出函数

建议新增一个独立实现文件，例如：

- `src/common/ber_window_probe.cpp`

职责：

- 聚合多 seed 的 per-window 统计
- 计算 `BER_agg[w]`
- 计算 `BER_ge_k`
- 计算 `delta(k)`
- 选出 `recommended_start_window`
- 导出：
  - `per_seed_per_window.csv`
  - `aggregated_window_ber.csv`
  - `start_window_candidates.csv`

这样能把“仿真运行”和“统计分析”分开。

### 16.2 建议修改的代码

#### 1. 修改 pipeline 结果结构，允许携带 window 级统计

建议修改：

- [pipeline_runner.hpp](/home/zsr71/projects/newcode/include/newcode/pipeline_runner.hpp)
- [pipeline_runner.cpp](/home/zsr71/projects/newcode/src/common/pipeline/pipeline_runner.cpp)

目标：

- 在不影响现有 BER 汇总输出的前提下，增加一条可选的 window 级统计结果链

建议做法：

- 在 `PipelineResult` 中新增一个 vector，例如：
  - `std::vector<WindowBerSample> window_ber_samples;`

注意：

- 这部分应当是“附加统计”
- 不要替换现有 `pre_fec / post_fec` 总体 BER 结果

#### 2. 修改 BER 统计代码，支持按 window 输出原始错误数

建议修改：

- [ber.cpp](/home/zsr71/projects/newcode/src/rx/ber/ber.cpp)

当前这里主要输出整体 BER。第一版可在不改变现有接口语义的前提下，新增一个辅助接口，例如：

- 保留原有 `compute_ber(...)`
- 新增：
  - `compute_ber_per_window(...)`
  - 或 `compute_ber_window_samples(...)`

职责：

- 仍使用和当前 BER 一致的 bit 对齐方式
- 只是把统计粒度从“整体”细化到“按 window”

这样可以保证：

- probe app 的 window 统计口径
- 与现有正式 BER 统计口径

保持一致。

#### 3. 修改 app 构建配置

建议修改：

- [CMakeLists.txt](/home/zsr71/projects/newcode/CMakeLists.txt)

目标：

- 把新的 `ofec_ber_window_probe` 编进工程

这是很小的改动，但第一版必须补上。

### 16.3 建议尽量不改的代码

第一版建议尽量不要修改这些部分：

- `apps/ofec_single.cpp`
- `apps/ofec_sweep.cpp`
- `apps/ofec_sweep3.cpp`
- 正式 runner 的默认 BER 裁剪逻辑
- warm-up 相关正式默认参数

原因是：

- 这次工作的目标是“先测、先观察、先定规则”
- 不是先改线上默认行为

### 16.4 推荐的实现顺序

建议按下面顺序落地：

1. 在 BER 层新增按 window 的统计接口
2. 在 `PipelineResult` 里挂上 window 级样本
3. 新增 `apps/ofec_ber_window_probe.cpp`
4. 实现原始 CSV 导出
5. 实现多 seed 聚合和 `BER_ge_k` 计算
6. 实现 `recommended_start_window` 自动推荐

这样做的好处是：

- 每一步都能单独验证
- 中间产物清楚
- 即使后面的推荐规则还没定死，前面的原始数据也已经可用了
