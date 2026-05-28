# OFEC Hybrid SISO Backfill 优先级策略设计

本文描述 hybrid prepass 中 `HYBRID_SISO_BACKFILL_MODE` 的一个新策略设想。

目标很明确：

- 当 `SISO budget` 宽松时，尽量把 `SISO` 用满。
- 当 `SISO budget` 紧张时，尽量把 `SISO` 留给更难的行。

本文只定义策略语义、变量含义和预期行为，不涉及代码修改细节。

## 1. 当前逻辑回顾

当前代码里，`HYBRID_SISO_BACKFILL_MODE` 的核心逻辑位于：

- `src/rx/ofec/detail/ofec_tile_dispatch.ipp`

当前已落地模式有两个：

- `Disabled`
- `TwoErrorOnly`

其中：

- `Disabled`
  - hybrid classifier 一旦返回 `hard_ok=true`，该行立即转成 `HardFinish`
- `TwoErrorOnly`
  - 只有 `HybridRowClass::TwoMain` 会先进入 `deferred_candidates`
  - tile 扫描结束后，再根据 `soft_rows_before_deferred_accept` 和 `siso_budget` 决定要不要把其中一部分收回成 `HardFinish`

当前回收公式可以概括为：

```text
max_deferred_accept =
    max(soft_rows_before_deferred_accept - siso_budget, 0)

deferred_accept_count =
    min(max_deferred_accept, deferred_candidates.size())
```

其中：

- `soft_rows_before_deferred_accept`
  - 当前所有 `tag == SoftDecode` 的行数
  - 既包含 deferred 的 `TwoMain`
  - 也包含 hard-finish 失败后继续留在 soft path 的行
- `siso_budget`
  - 当前 tile 的 soft decode 预算

所以当前 `TwoErrorOnly` 的语义是：

- 预算不紧时，先把 `TwoMain` 留在 soft path，尽量把 `SISO` 用满
- 预算超出时，优先从 deferred 的 `TwoMain` 里收回一部分

## 2. 新策略动机

当前 `TwoErrorOnly` 的优点是简单，但它只有一种“可回收 hard-finish”：

- `TwoMain`

这意味着：

- `OneMain` / `OneMainPlusParity` 即使已经能稳定 hard-finish，也会立即退出 soft path
- 当 `SISO budget` 宽松时，系统无法用这些“容易行”去补满 soft path
- 当 `SISO budget` 紧张时，系统也无法表达“先踢一错，再踢二错”的资源偏好

如果目标是：

- 预算宽松时尽量用满 `SISO`
- 预算紧张时把 `SISO` 留给更难的行

那么更合适的偏好顺序应当是：

- `OneMain`
- `OneMainPlusParity`
- `TwoMain`
- 其他 hard-finish 失败行

这里的含义不是“这些类一定更重要”，而是：

- 一错行已经足够容易，最适合优先收回成 `HardFinish`
- 二错行比一错行更值得保留在 soft path 里
- hard-finish 失败行默认视为更难行，优先保留给 `SISO`

## 3. 新模式定义

建议新增一个 backfill 模式，本文暂称：

- `OneAndTwoErrorPriority`

该模式的语义如下：

- deferred 候选扩展为：
  - `OneMain`
  - `OneMainPlusParity`
  - `TwoMain`
- 同类内部保持当前逐行扫描顺序
- tile 扫描结束后，若 `soft_rows_before_deferred_accept <= siso_budget`
  - 不回收任何 deferred 候选
  - 尽量把 `SISO` 用满
- 若 `soft_rows_before_deferred_accept > siso_budget`
  - 先回收 `OneMain`
  - 再回收 `OneMainPlusParity`
  - 最后回收 `TwoMain`

这一定义表达的是一个明确的资源优先级：

```text
最先踢出：
  OneMain

其次踢出：
  OneMainPlusParity

最后踢出：
  TwoMain
```

也可以把它合并理解为：

```text
先踢所有一错类
再踢二错类
同类内部按扫描顺序
```

## 4. 变量语义

本文延续当前代码里的两个关键量：

### 4.1 `soft_rows_before_deferred_accept`

它的定义不变：

- 统计当前所有 `tag == SoftDecode` 的行
- 不是“只统计 hard-fail”
- 也不是“只统计 deferred 的二错”

因此在新策略下，它会同时包含：

- deferred 的 `OneMain`
- deferred 的 `OneMainPlusParity`
- deferred 的 `TwoMain`
- 其他 hard-finish 失败后保留在 soft path 的行

### 4.2 `siso_budget`

它的定义也不变：

- 当前 tile 最多允许多少行进入 soft decode / Chase

## 5. 回收规则

新模式建议保留当前“只在 soft 超预算时才回收”的总体结构，只改变回收候选集合和优先级。

### 5.1 计算超额数

```text
overflow = max(soft_rows_before_deferred_accept - siso_budget, 0)
```

含义：

- `overflow == 0`
  - soft 行数没有超过预算
  - 一个 deferred 都不回收
- `overflow > 0`
  - soft 行数超预算
  - 需要从 deferred 池里收回最多 `overflow` 个

### 5.2 分级回收

建议按下面顺序回收：

```text
step 1:
  从 deferred OneMain 中按扫描顺序回收

step 2:
  若仍有 overflow
  从 deferred OneMainPlusParity 中按扫描顺序回收

step 3:
  若仍有 overflow
  从 deferred TwoMain 中按扫描顺序回收
```

### 5.3 同类内部顺序

同类内部不新增排序规则，保持当前实现风格：

- 谁先在 tile 行扫描时进入 deferred
- 谁先被回收

这样有两个好处：

- 第一版行为最容易解释
- 不需要引入新的 reliability score 或额外 side data

## 6. 数字例子

下面给几个具体例子。

### 6.1 预算宽松：不回收任何 deferred

假设：

- `siso_budget = 16`
- `soft_rows_before_deferred_accept = 14`
- 这 14 行里包含：
  - `3` 个 `OneMain`
  - `2` 个 `OneMainPlusParity`
  - `4` 个 `TwoMain`
  - `5` 个 hard-finish 失败行

则：

```text
overflow = max(14 - 16, 0) = 0
```

结果：

- 一个都不回收
- 这 14 行全部继续留在 soft path

### 6.2 预算刚好：仍然不回收

假设：

- `siso_budget = 16`
- `soft_rows_before_deferred_accept = 16`

则：

```text
overflow = max(16 - 16, 0) = 0
```

结果：

- 一个都不回收
- soft path 正好打满预算

### 6.3 预算紧张：先回收一错

假设：

- `siso_budget = 16`
- `soft_rows_before_deferred_accept = 19`
- 当前 soft 行构成为：
  - `2` 个 `OneMain`
  - `1` 个 `OneMainPlusParity`
  - `4` 个 `TwoMain`
  - `12` 个 hard-finish 失败行

则：

```text
overflow = 19 - 16 = 3
```

回收顺序：

1. 先回收 `2` 个 `OneMain`
2. 还差 `1` 个，再回收 `1` 个 `OneMainPlusParity`
3. `TwoMain` 不动

最终 soft path 保留：

- `0` 个 `OneMain`
- `0` 个 `OneMainPlusParity`
- `4` 个 `TwoMain`
- `12` 个 hard-finish 失败行

总数正好 `16`

### 6.4 一错不够时，再回收二错

假设：

- `siso_budget = 16`
- `soft_rows_before_deferred_accept = 19`
- 当前 soft 行构成为：
  - `1` 个 `OneMain`
  - `0` 个 `OneMainPlusParity`
  - `4` 个 `TwoMain`
  - `14` 个 hard-finish 失败行

则：

```text
overflow = 19 - 16 = 3
```

回收顺序：

1. 先回收 `1` 个 `OneMain`
2. 还差 `2` 个，再从 `TwoMain` 中按扫描顺序回收 `2` 个

最终 soft path 保留：

- `0` 个 `OneMain`
- `0` 个 `OneMainPlusParity`
- `2` 个 `TwoMain`
- `14` 个 hard-finish 失败行

总数正好 `16`

## 7. 预期收益

这个策略的主要收益是：

- 预算宽松时，`OneMain` / `OneMainPlusParity` / `TwoMain` 都可以先留在 soft path，帮助吃满 `SISO`
- 预算紧张时，可以显式表达“先踢容易行，再保留难行”的优先级
- 第一版不需要引入新的置信度排序，行为仍然容易解释

它体现的是一种工程偏好：

- 已经容易 hard-finish 的行，soft decode 的边际价值更低
- 越难的行，越值得优先保留 `SISO`

## 8. 风险和边界

需要明确一点：

- “hard-finish 失败行”不等于“最值得 soft 的行”
- 它们只是“当前 classifier / hard prepass 没有直接收掉的行”

因此，新策略本质上是在用“分类难度”近似替代“soft decode 预期收益”。

这通常是合理的，但不保证始终最优。可能存在以下情况：

- 某些 hard-finish 失败行其实很差，继续占 `SISO` 收益不高
- 某些 `TwoMain` 行即使 hard-finish 成功，若继续跑 soft，也可能得到更好的外信息

所以这套策略更适合作为：

- 第一版优先级规则
- 一个清晰、稳定、容易验证的工程基线

而不是最终最优策略。

## 9. 建议的模式语义

如果后续要把它落到枚举层，建议语义写得非常明确：

- `Disabled`
  - 不做 backfill
- `TwoErrorOnly`
  - 只 deferred `TwoMain`
  - 超预算时只回收 `TwoMain`
- `OneAndTwoErrorPriority`
  - deferred `OneMain`、`OneMainPlusParity`、`TwoMain`
  - 超预算时按
    - `OneMain`
    - `OneMainPlusParity`
    - `TwoMain`
    的顺序回收
  - 同类内部按扫描顺序

## 10. 伪代码

下面给出一版与当前 `run_hybrid_prepass()` 结构对齐的伪代码。

```text
deferred_one_main = []
deferred_one_main_plus_parity = []
deferred_two_main = []

for each soft row:
    run hybrid classifier

    if hard_ok == false:
        keep row as SoftDecode
        continue

    if mode == Disabled:
        apply_hard_finish_to_plan(row)
        continue

    if mode == TwoErrorOnly:
        if hard_class == TwoMain:
            deferred_two_main.push(row)
            continue
        apply_hard_finish_to_plan(row)
        continue

    if mode == OneAndTwoErrorPriority:
        if hard_class == OneMain:
            deferred_one_main.push(row)
            continue
        if hard_class == OneMainPlusParity:
            deferred_one_main_plus_parity.push(row)
            continue
        if hard_class == TwoMain:
            deferred_two_main.push(row)
            continue
        apply_hard_finish_to_plan(row)
        continue

overflow = max(soft_rows_before_deferred_accept - siso_budget, 0)

reclaim from deferred_one_main in scan order
reclaim from deferred_one_main_plus_parity in scan order
reclaim from deferred_two_main in scan order

rebuild_soft_candidate_rows()
```

## 11. 结论

如果目标是：

- 预算宽松时尽量用满 `SISO`
- 预算紧张时优先把 `SISO` 留给更难行

那么把 backfill 候选从“只含 `TwoMain`”扩展为“`OneMain` + `OneMainPlusParity` + `TwoMain`”，并在超预算时按“一错优先回收、二错后回收”的顺序处理，是一个合理且实现风险较低的第一版策略。
