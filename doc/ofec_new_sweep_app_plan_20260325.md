# 新增仿 `apps/ofec_sweep.cpp` 的 App 实现方案

## 目标

在 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 之外，再新增一个独立的 app 入口，用于承载一套新的 sweep 配置，但尽量复用现有的 sweep 基础设施：

- 复用 [`include/newcode/ofec_sweep_runner.hpp`](/home/zsr71/projects/newcode/include/newcode/ofec_sweep_runner.hpp)
- 复用 [`src/ofec_sweep/`](/home/zsr71/projects/newcode/src/ofec_sweep)
- 复用现有的 `SweepParameterConfig -> run_sweep(config)` 执行链路

这类新 app 的定位不是再发明一套新的 runner，而是：

- 提供一个新的“参数装配入口”
- 用一组和 `ofec_sweep.cpp` 不同的默认配置组织实验
- 保持构建、运行、输出、日志口径与现有 `ofec_sweep` 兼容

## 适用场景

适合新增这个 app 的情况：

- 想保留 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 现有用途，不希望直接改乱主入口
- 想做一套独立的参数搜索实验，例如只扫部分轴，或固定某些参数
- 想用另一组默认 `alpha/beta`、`early-stop`、`MUX`、`seed`、`Eb/N0` 配置
- 想让实验入口更聚焦，避免一个 app 顶部堆太多互相无关的常量

不适合新增 app 的情况：

- 只是改一两个常量做一次性实验
- 只是想在现有 `ofec_sweep.cpp` 里多加一个候选值
- 需要新的执行语义而不是新的参数入口

如果需求只是“换一组默认参数”，原则上优先复用现有 runner，新 app 只负责装配配置。

## 建议文件

建议新增：

- `apps/<new_app_name>.cpp`

可选新增：

- 在 [`CMakeLists.txt`](/home/zsr71/projects/newcode/CMakeLists.txt) 里注册新的可执行目标
- 在 [`doc/`](/home/zsr71/projects/newcode/doc) 下补充该 app 的用途说明

其中 `<new_app_name>` 建议使用能体现用途的名字，例如：

- `ofec_sweep_es.cpp`
- `ofec_sweep_fixed_decoder.cpp`
- `ofec_sweep_mux_reconfig.cpp`
- `ofec_sweep_ablation.cpp`

避免使用没有语义的信息名，例如 `ofec_sweep3.cpp`，除非确实只是临时实验入口。

## 设计原则

### 1. 新 app 只做“配置装配”

新 app 应尽量像 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 一样：

- 顶部定义一组用户可调常量
- `main()` 中把这些常量组装成 `ofec_sweep::SweepParameterConfig`
- 最后直接调用 `ofec_sweep::run_sweep(config)`

不要在 app 里重复实现：

- scenario 展开
- 多线程执行
- CSV 输出
- 最优结果汇总
- 参数校验

这些都已经在 runner 层具备。

### 2. 和 `ofec_sweep.cpp` 保持同一套参数语义

新 app 中的字段命名、注释口径、默认值组织方式，建议与 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 保持一致。这样后续迁移配置时成本最低，也更容易横向比较。

建议继续保留以下分区：

- 发端参数
- 信道参数
- 量化参数
- 解码参数
- early-stop 参数
- MUX 参数
- debug 参数

### 3. 差异只体现在“默认配置”和“暴露的扫描轴”

新 app 与 `ofec_sweep.cpp` 的差异最好体现在：

- 只暴露本实验真正关心的扫描维度
- 其他无关参数固定
- 显式说明哪些参数被固定，哪些参数会展开扫描

这样能避免出现“看起来是新 app，实际上和老 app 只是复制粘贴了一份但更难维护”的问题。

## 推荐实现结构

建议直接以 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 为模板，保留它的总体结构：

```cpp
#include <vector>

#include "mux_bypass_edges.hpp"
#include "newcode/ofec_sweep_runner.hpp"
#include "newcode/utils/linspace.hpp"

namespace {

// 1. 顶部常量区
// 2. 按“发端 / 信道 / 量化 / 解码 / early-stop / MUX / debug”分段

}  // namespace

int main() {
  const auto& selected_mux_bypass_edges =
      app_mux::bypass_edges_for_scheme(kMuxBypassScheme);

  ofec_sweep::SweepParameterConfig config;

  // 1. 基础链路配置
  // 2. 早停固定配置
  // 3. 扫描候选
  // 4. MUX / 调度配置
  // 5. seed 配置
  // 6. Eb/N0 配置
  // 7. debug trace 配置
  // 8. base_params 基础字段

  return ofec_sweep::run_sweep(config);
}
```

推荐做法是“先复制结构，再删掉不用的轴”，而不是从空文件重新拼。

## 新 app 应明确回答的 5 个问题

在开始写代码前，先把这 5 个问题定下来：

1. 新 app 的实验目标是什么？
2. 它和 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 的核心差异是什么？
3. 哪些参数要固定，哪些参数要扫描？
4. 输出结果里最关心的是 BER、FER、early-stop 比例，还是最优参数组合？
5. 是否需要显式 `explicit_patterns`，还是继续使用 `start/step` 方式生成 `ALPHA_LIST / beta_list`？

如果这 5 个问题不先定清楚，最后很容易写成“又一个什么都能扫，但没人知道为什么存在的 app”。

## 参数层面的建议

### 1. `decoder_name` 和候选列表

如果新 app 只服务一种 decoder，建议：

- 固定 `config.decoder_name`
- `config.decoder_name_candidates` 置空

如果新 app 的目标是比较 decoder 变体，才开启 `decoder_name_candidates`。

### 2. `explicit_patterns` 和 `alpha/beta` 网格

如果实验重点是“比较一批手工设计的 `ALPHA_LIST / beta_list`”：

- 优先使用 `config.explicit_patterns`

如果实验重点是“从简单规则生成一批参数”：

- 使用 `alpha_start_candidates`
- 使用 `alpha_step_candidates`
- 使用 `beta_start_candidates`
- 使用 `beta_step_candidates`

不要两套都大量展开，否则 scenario 数量会膨胀得很快。

### 3. early-stop 相关轴

如果新 app 主要研究 early-stop，建议把扫描轴集中在：

- `early_stop_condition_candidates`
- `early_stop_action_candidates`
- `early_stop_v2_llr_abs_threshold_candidates`
- `early_stop_v2_max_unreliable_bits_candidates`
- `early_stop_action_beta_start_candidates`
- `early_stop_action_beta_step_candidates`

同时尽量固定 Chase 和 MUX 参数，避免结论混杂。

### 4. MUX 相关轴

如果新 app 主要研究调度策略，建议聚焦：

- `siso_active_list`
- `mux_group_g`
- `mux_enable_reconfig`
- `mux_bypass_scheme`

同时把 `alpha/beta` 和 early-stop 条件固定住。

## 实现步骤

### 第一步：复制入口骨架

从 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 复制出一个新文件：

- 保留 include
- 保留 `namespace` 顶部常量区结构
- 保留 `main()` 中 `SweepParameterConfig` 的装配顺序

### 第二步：删掉与新目标无关的扫描轴

例如如果新 app 只研究 early-stop，就删除或固定：

- `decoder_name_candidates`
- `chase_topk_keep_candidates`
- `chase_group_minima_bits_candidates`
- `explicit_patterns`

反过来，如果新 app 只研究 `alpha/beta`，就固定：

- `early_stop_condition_mode`
- `early_stop_action_mode`
- `early_stop_*_candidates`

### 第三步：补充新的默认常量

把该实验真正关心的配置显式写在顶部常量区，例如：

- 新的 `Eb/N0` 范围
- 新的 seed 口径
- 新的 `SISO_ACTIVE_LIST`
- 新的 `kExplicitAlphaBetaSets`
- 新的 early-stop 配置

### 第四步：确认 `base_params` 与顶层配置一致

需要特别留意这些字段不要漏：

- `config.base_params.BITGEN_RANDOM_BITS`
- `config.base_params.NORMALIZE_KNOWN_PREFIX_TAIL`
- `config.base_params.ENABLE_EARLY_STOP`
- `config.base_params.LLR_CLIP_RATIO`
- `config.base_params.LLR_BITS`
- `config.base_params.SISO_ACTIVE_LIST`

如果顶层 `config.xxx` 和 `base_params.xxx` 语义重复，必须保持一致，避免运行时出现“看起来设置了，但实际没生效”的误解。

### 第五步：接入构建

在 [`CMakeLists.txt`](/home/zsr71/projects/newcode/CMakeLists.txt) 中新增可执行目标，使其和现有 app 一样可以单独编译和运行。

### 第六步：做最小验证

至少验证：

1. 能成功编译
2. 能成功启动并跑完一个最小 sweep
3. 输出 CSV 与控制台结果正常
4. scenario 数量符合预期
5. 关键参数确实出现在结果命名或日志中

## 最小验收标准

满足以下条件即可认为新 app 落地完成：

1. `cmake --build build --target <new_app_name> -j` 可以通过
2. `./build/<new_app_name>` 可以跑通
3. 会调用 [`ofec_sweep::run_sweep()`](/home/zsr71/projects/newcode/include/newcode/ofec_sweep_runner.hpp)
4. 不新增新的 sweep runner 分支
5. 新 app 的差异点清晰，且不是 `ofec_sweep.cpp` 的无意义复制

## 推荐的文档补充

新 app 建好后，建议再补两类文档：

- 运行说明：如何编译、如何运行、输出文件在哪
- 设计说明：它和 [`apps/ofec_sweep.cpp`](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp) 的区别

否则几周后再回头看，很难知道“为什么当时需要再开一个 app”。

## 一个务实的落地建议

如果现在需求还不够具体，建议先按下面的方式推进：

1. 先确定新 app 名字
2. 先明确它只研究哪一类问题
3. 先把扫描轴压到 1 到 3 个
4. 先做一个能跑通的小入口
5. 跑通后再逐步放开更多候选轴

这样比一开始就做成“第二个万能 `ofec_sweep.cpp`”更稳，也更容易维护。
