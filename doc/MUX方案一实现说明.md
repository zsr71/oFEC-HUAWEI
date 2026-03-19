# MUX方案一实现说明

## 1. 目标

方案一的目标是实现这样一种行为：

- 当某个 tile 内可用的 SISO 数量 **不小于** 当前需要处理的 code 数量时，不再做 MUX 调度
- 当某个 tile 内可用的 SISO 数量 **小于** 当前需要处理的 code 数量时，才进入现有的 MUX 分组 / 重配置逻辑

也就是说，这个方案想表达的是：

> 只有在“资源不够”的情况下，MUX 才真正介入；资源够的时候，MUX 应该退化为空操作

这和当前的实验直觉是一致的：

- `32 code -> 32 SISO` 时，不需要调度
- `32 code -> 16 SISO` 时，才需要调度、分组、旁路

---

## 2. 为什么这是最小改动方案

方案一不改变现有参数结构，也不引入新的按 tile 配置项。

它保留当前所有顶层参数：

- `SISO_ACTIVE_LIST`
- `MUX_GROUP_G`
- `MUX_ENABLE_RECONFIG`
- `MUX_EXTRA_BYPASS_EDGES`

只是在 tile 内真正执行 MUX 之前，增加一道判断：

- 如果当前 tile 的 SISO 预算已经够用，就直接跳过 MUX 调度

因此它的特点是：

- 改动范围小
- 风险低
- 不会破坏当前 sweep 和 single 的顶层接口
- 非常适合先验证“只有资源不足时才需要 MUX”这个想法

---

## 3. 现在的流程是什么

当前 tile 内的大致流程是：

1. 先根据 early-stop 结果构造 `mux_state`
2. 不管资源是否足够，都会继续进入：
   - 普通 grouped budget 路径，或
   - reconfig 路径
3. 最后得到哪些 code 本轮允许跑 SISO，哪些会变成 `Unscheduled`

也就是说，当前实现里：

- `MUX_GROUP_G`
- `MUX_ENABLE_RECONFIG`
- `MUX_EXTRA_BYPASS_EDGES`

这些逻辑只要进了 tile，就会参与判断。

这就导致一个现象：

- 即使某个 tile 的 `siso_active_for_tile >= rows_to_decode`
- 理论上根本没有资源冲突
- 程序仍然会走 MUX 的一套流程

从实验角度看，这不是最理想的行为。

---

## 4. 方案一的核心改法

方案一的核心就是在 tile 级增加一个“资源是否紧张”的判断。

推荐判断条件：

- `siso_active_for_tile >= rows_to_decode`

这里：

- `rows_to_decode` 表示这个 tile 当前真正参与 Chase / 行解码的 code 数
- `siso_active_for_tile` 表示这个 tile 当前可用的 SISO 数量

因此：

- 如果 `siso_active_for_tile >= rows_to_decode`
  - 说明资源够用
  - 所有需要解码的 code 都可以直接保留
  - MUX 不需要做任何裁剪或重配置

- 如果 `siso_active_for_tile < rows_to_decode`
  - 说明资源不足
  - 这时才进入现有的 grouped / reconfig MUX 逻辑

---

## 5. 建议改动位置

### 5.1 主要改动点：`ofec_tile_impl.ipp`

最核心的改动位置是：

- [src/rx/ofec/detail/ofec_tile_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/detail/ofec_tile_impl.ipp)

原因是：

- 这里已经拿到了 `rows_to_decode`
- 这里也已经拿到了 `siso_active_for_tile`
- 并且当前的 `mux_state` 正是在这里被构造和调度的

所以最自然的实现方式是：

1. 先完成 `early_stop_stats`
2. 先构造 `mux_state`
3. 再判断当前 tile 是否“资源足够”
4. 如果资源足够：
   - 直接跳过后面的 grouped / reconfig 调度
5. 如果资源不足：
   - 继续走现有 MUX 流程

这意味着当前代码里这两段逻辑：

- `schedule_scheme_c_staged_cpp(...)`
- `apply_siso_budget_grouped(...)`

都应该只在“资源不足”分支下执行。

---

### 5.2 可选辅助改动：增加一个小的辅助判断函数

为了让代码更清晰，可以把“资源是否足够”的判断单独抽成一个小函数，例如：

- `mux_is_effectively_needed(rows_to_decode, siso_active_for_tile)`

这个函数可以放在：

- `ofec_tile_impl.ipp` 内部匿名区域

也可以单独放到一个小的 MUX helper 文件里。

它的作用只是让主流程更清楚，例如：

- `true` 表示当前 tile 需要 MUX
- `false` 表示当前 tile 不需要 MUX

这一步不是必须的，但会让后续继续扩展时更整洁。

---

## 6. 具体行为应该是什么

方案一落地后，建议行为定义成下面这样。

### 6.1 资源足够时

条件：

- `siso_active_for_tile >= rows_to_decode`

行为：

- 保持 `mux_state` 中原本的 `NeedSiso` / `EarlyStopped` 标记
- 不再额外把任何 `NeedSiso` 行改成 `Unscheduled`
- 不走 grouped budget
- 不走 reconfig

这意味着：

- 所有本该做 SISO 的行都继续保留
- MUX 在这个 tile 上退化为空操作

---

### 6.2 资源不足时

条件：

- `siso_active_for_tile < rows_to_decode`

行为：

- 继续沿用当前已有逻辑
- 如果 `MUX_ENABLE_RECONFIG=true`
  - 走 staged reconfig 调度
- 如果 `MUX_ENABLE_RECONFIG=false`
  - 走 grouped budget 裁剪

这意味着：

- 方案一不会重写已有 MUX 算法
- 只是增加一个“是否需要进入 MUX”的前置判断

---

## 7. 对现有参数的影响

方案一的优点之一是：**不需要新增顶层参数**。

现有参数的含义保持不变：

- `SISO_ACTIVE_LIST`
  - 继续表示每个 tile 的 SISO 预算

- `MUX_GROUP_G`
  - 继续表示分组粒度

- `MUX_ENABLE_RECONFIG`
  - 继续表示是否启用重配置版调度

- `MUX_EXTRA_BYPASS_EDGES`
  - 继续只在 reconfig 路径下生效

唯一变化是这些参数的生效前提变成：

> 只有当当前 tile 的 SISO 预算不足时，这些 MUX 参数才真正发挥作用

这其实更符合直觉。

---

## 8. 对实验的意义

如果采用方案一，你就能自然得到下面这种实验效果：

- `SISO_ACTIVE_LIST = {32,32,32,32}`
  - 所有 tile 都不触发有效 MUX

- `SISO_ACTIVE_LIST = {32,32,32,16}`
  - 前三个 tile 不触发有效 MUX
  - 最后一个 tile 才触发有效 MUX

- `SISO_ACTIVE_LIST = {16,16,16,16}`
  - 所有 tile 都触发有效 MUX

这比当前“所有 tile 都统一走 MUX 流程”更符合你想表达的实验语义。

---

## 9. 优点和限制

### 优点

- 改动最小
- 不需要改顶层参数接口
- 不影响已有 `ofec_single` / `ofec_sweep` 配置结构
- 非常适合快速验证“资源不足时才需要 MUX”这个想法

### 限制

- `MUX_GROUP_G` 仍然是全 tile 共用的
- `MUX_ENABLE_RECONFIG` 仍然是全 tile 共用的
- `MUX_BYPASS_SCHEME` 仍然是全 tile 共用的

也就是说，方案一虽然能解决“资源足够时不该调度”的问题，但还不能实现：

- 前 3 个 tile 不分组
- 最后 1 个 tile 才分组

如果后续你想做到这种更细粒度的控制，就需要进入方案二。

---

## 10. 推荐的实现顺序

建议按下面顺序做。

1. 先在 `ofec_tile_impl.ipp` 增加“资源是否足够”的判断
2. 确认资源足够时，MUX 完全跳过
3. 用几组典型配置验证行为：
   - `{32,32,32,32}`
   - `{32,32,32,16}`
   - `{16,16,16,16}`
4. 确认日志和 BER 结果符合预期后，再决定是否继续做方案二

---

## 11. 一句话总结

方案一的本质是：

> 不改变现在的 MUX 算法，只增加一个前置判断，让 MUX 只在“SISO 不够用”的 tile 上真正生效

这是最小改动、最适合先验证实验想法的一种实现方式。

