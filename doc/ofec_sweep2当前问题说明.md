# `ofec_sweep2` 当前问题说明

## 1. 文档目的

这份说明用于记录 `apps/ofec_sweep2.cpp` 当前实现里已经看到的主要问题。  
重点不是代码细节，而是说明：

- 现在有哪些现象
- 这些现象为什么会出现
- 它们会造成什么影响
- 后续适合怎么修

---

## 2. 当前定位

`ofec_sweep2` 不是直接复用标准 `ofec_sweep` 的场景生成逻辑，而是：

- 自己构造参数形状
- 自己生成 `ExplicitAlphaBetaPattern`
- 自己手工拼 `SweepScenario`
- 再调用公共的并行运行接口

这种写法的好处是灵活，适合做“两阶段筛选”类实验。  
但代价是：它没有完整复用 `ofec_sweep` 现在的参数灌入逻辑，因此很容易出现“真实运行参数”和“记录到日志/CSV 的参数”不一致的问题。

---

## 3. 已看到的主要问题

### 3.1 `early_stop_action_sign_beta_list` 没有填，导致场景被跳过

现在 `ofec_sweep2` 里生成 pattern 时，只填了：

- `alpha_list`
- `beta_list`

但没有填：

- `early_stop_action_sign_beta_list`

而公共 runner 在提交场景前，会检查下面三条列表长度是否都等于 `TILES_PER_WIN`：

- `alpha_list`
- `beta_list`
- `early_stop_action_sign_beta_list`

当前 `ofec_sweep2` 里 `kTilesPerWindow = 2`，因此它要求三条列表长度都必须是 `2`。  
如果 `early_stop_action_sign_beta_list` 为空，就会出现这种告警：

`[WARN] Scenario '...' skipped due to list size mismatch (expected 2)`

### 影响

- 场景会被直接跳过
- 后续阶段可用候选数减少
- 如果所有场景都被跳过，阶段运行会直接失败

### 根本原因

`ofec_sweep2` 构造的 pattern 信息不完整，没有满足当前公共 runner 对 tile 级列表的长度要求。

---

### 3.2 CSV 中记录出的部分参数，并不一定是这次真正运行的参数

`ofec_sweep2` 写 CSV 时，调用的是公共的 `Extended CSV` 写行逻辑。  
这套 CSV 会写很多字段，例如：

- `chase_group_minima_bits`
- early-stop 条件/动作相关参数
- `early_stop_action_hard_llr_mag`

但 `ofec_sweep2` 手工构造 `SweepScenario` 时，并没有把这些字段按当前配置完整写进去。  
于是这些字段就会落回结构体默认值。

### 影响

- CSV 里看起来“每列都有值”
- 但其中一部分值未必是这次实验真正采用的值
- 后续分析时，容易把某些结果错误归因到这些字段上

### 根本原因

`ofec_sweep2` 没有走标准 `build_scenarios()` 路径，而是手工构造 scenario。  
一旦手工构造时漏填字段，日志/CSV 就会和实际运行状态脱节。

---

### 3.3 `chase_n_test` 现在是按 `2^L` 手工推出来的，不是显式取真实配置

`ofec_sweep2` 里构造 scenario 时，当前是按：

- `chase_L = base_params.CHASE_L`
- `chase_n_test = 1 << chase_L`

来写入 `SweepScenario`。

这意味着它默认认为：

- `CHASE_NTEST` 永远等于 `2^CHASE_L`

但这只适用于“满枚举 test pattern”的情况。  
如果以后想单独改 `CHASE_NTEST`，或者做“只取部分候选”的试验，这里就会和真实配置不一致。

### 影响

- CSV 和日志里记录的 `chase_n_test` 可能是错的
- 实际运行的是一个值，记录出来的是另一个值

### 根本原因

`ofec_sweep2` 手工推导了一个“推测值”，而不是直接从最终运行参数里取值。

---

### 3.4 `SweepParameterConfig` 里的很多字段在 `ofec_sweep2` 中并不会自动生效

`ofec_sweep2` 仍然构造了完整的 `SweepParameterConfig`，但并没有把它完整交给标准 `ofec_sweep` 的场景展开逻辑。  
它只是把其中一部分字段抽出来，自己再拼成 scenario。

因此会出现这种情况：

- 顶层配置里看起来有很多字段
- 但不是所有字段都会自动进入最终 scenario

### 影响

容易出现误解，例如：

- 改了某个 `config` 字段
- 以为它会像 `ofec_sweep` 那样自动进入所有场景
- 实际上并没有

### 根本原因

`ofec_sweep2` 现在本质上是“半复用 `ofec_sweep` 基础设施”的写法，不是完整复用。

---

## 4. 这些问题的共同本质

`ofec_sweep2` 当前的主要风险，不是“程序完全跑不起来”，而是：

- 有些场景会因为列表长度不完整被跳过
- 有些结果文件里记录出来的参数，不一定就是这次真实跑的参数

也就是说，它最大的风险是：

**结果解释容易出错。**

这比单纯的崩溃更麻烦，因为实验可能能跑完，但后面做分析时会被错误元数据带偏。

---

## 5. 建议的修复方向

### 方案 A：最小修复

目标是先让 `ofec_sweep2` 的当前实验能稳定跑、结果也不至于误导。

建议做法：

1. 在 pattern 里把 `early_stop_action_sign_beta_list` 补齐到 `TILES_PER_WIN`
2. 在手工构造 `SweepScenario` 时，把当前真正会写进 CSV 的字段都补全
3. `chase_n_test` 不再手工推导，直接取最终实际配置值

### 优点

- 改动小
- 适合快速恢复实验可用性

### 限制

- `ofec_sweep2` 仍然是手工拼 scenario
- 以后如果公共 `SweepScenario` 再加字段，还要继续人工同步

---

### 方案 B：让 `ofec_sweep2` 更多复用标准场景展开

目标是减少“手工拼 scenario”带来的字段漏填问题。

思路是：

- `ofec_sweep2` 主要保留“两阶段筛选”的实验逻辑
- 但具体的 scenario 元数据生成尽量复用标准 `build_scenarios()` 或同一套填充逻辑

### 优点

- 和 `ofec_sweep` 的行为更一致
- 以后新增参数时，不容易漏接

### 限制

- 改动比方案 A 大
- 需要重新整理 `ofec_sweep2` 现在的实验组织方式

---

## 6. 当前最值得优先修的点

如果只按实验优先级排序，我建议先修下面三个：

1. `early_stop_action_sign_beta_list` 长度问题  
   否则场景会被直接跳过

2. `SweepScenario` 字段补全问题  
   否则 CSV 会记录错参数

3. `chase_n_test` 不要再用 `2^L` 手工推  
   直接用真实配置值

---

## 7. 总结

`ofec_sweep2` 当前不是“完全不可用”，而是处在一个比较典型的状态：

- 运行主流程能工作
- 但实验元数据和公共 runner 的约束没有完全对齐

因此它现在最需要做的，不是增加更多功能，而是先把：

- 场景输入完整性
- 结果记录一致性

这两件事收紧。

只有把这两个基础问题收好，后续两阶段筛选的结果才适合拿来做稳定比较。
