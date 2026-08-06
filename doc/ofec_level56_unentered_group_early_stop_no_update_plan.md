# Level 5/6 EarlyStop 组更新三模式实现方案

## 1. 目的

本文描述 Level 5/6 `Group4LoadSortedMultiround` 调度下，EarlyStop 行是否更新外信息的三种可切换模式，以及最终代码实现。

本功能只通过 `apps/ofec_sweep3.cpp` 中的一个常量切换模式。一次 sweep 只运行所选模式，不增加 scenario 维度，也不自动展开三份实验。

## 2. 调度背景

Level 5 和 Level 6 共 64 个 code，固定分为 16 组，每组 4 个：

```text
Group 0  = shared code 0-3
Group 1  = shared code 4-7
...
Group 15 = shared code 60-63
```

所有组累计最多使用 8 次 group entry。普通调度每选中一次组，该组最多规划一个 SISO code 和一个 HISO code。若组内仍有未调度的普通候选，该组可以在后续轮次再次进入。

EarlyStop code 不需要 HISO/SISO core，但是否允许它更新外信息，取决于本功能的模式。

## 3. 参数定义

参数采用枚举，不再使用原来的布尔开关：

```cpp
enum class Level56EarlyStopGroupUpdateMode : uint8_t {
  AllGroups = 0,
  EnteredGroupsOnly = 1,
  FillIdleEntries = 2
};
```

`Params` 中的默认值为：

```cpp
Level56EarlyStopGroupUpdateMode LEVEL56_EARLY_STOP_GROUP_UPDATE_MODE =
    Level56EarlyStopGroupUpdateMode::AllGroups;
```

默认使用 `AllGroups`，保证未显式配置该参数的程序保持历史行为。

## 4. 三种模式

### 4.1 `AllGroups`

所有 EarlyStop 命中行均执行 `EarlyStopAction`，不考虑所在组是否获得 group entry。

这是历史行为。

### 4.2 `EnteredGroupsOnly`

先按原调度规则完成普通 HISO/SISO 调度。只有至少获得过一次普通 group entry 的组，才允许组内 EarlyStop 行执行 `EarlyStopAction`。

没有进入过的组，其 EarlyStop 行改为：

```text
final_action = Unscheduled
produced = false
```

因此这些行不计算新的 EarlyStop 外信息，也不写回。

### 4.3 `FillIdleEntries`

先完整执行与原来相同的普通多轮 HISO/SISO 调度。普通候选全部处理完后，如果 8 次 group entry 尚未用满，则用剩余次数访问此前从未进入过的组。

补位选择规则固定为：

```text
按 group index 从小到大选择尚未进入的组
```

补位访问只允许该组的 EarlyStop 行更新外信息，不规划 HISO/SISO code，也不分配 HISO/SISO core。

普通候选始终优先。只要仍有普通 HISO/SISO 候选，剩余 entry 就继续用于原有多轮调度，不会提前用于 EarlyStop 补位。

## 5. `K` 不同取值下的行为

这里 `K` 表示第一轮开始前，包含普通非 EarlyStop 候选的非空组数量。

### 5.1 `K > 8`

第一轮已经用完 8 次 entry，没有补位空间：

```text
AllGroups:
  所有 16 组的 EarlyStop 行都更新。

EnteredGroupsOnly:
  只有普通调度选中的 8 组更新。

FillIdleEntries:
  与 EnteredGroupsOnly 相同，因为没有空闲 entry。
```

### 5.2 `K = 8`

第一轮恰好用完 8 次 entry：

```text
AllGroups:
  所有 16 组的 EarlyStop 行都更新。

EnteredGroupsOnly:
  第一轮进入的 8 组更新。

FillIdleEntries:
  与 EnteredGroupsOnly 相同，没有空闲 entry。
```

### 5.3 `0 < K < 8`

第一轮先让这 `K` 个非空组进入。之后原调度器继续处理这些组中剩余的普通候选，直至候选耗尽或累计用满 8 次 entry。

普通调度结束后：

```text
remaining = 8 - used_group_entries
```

三种模式分别为：

```text
AllGroups:
  不使用 remaining 做补位，但所有 16 组的 EarlyStop 行仍更新。

EnteredGroupsOnly:
  不使用 remaining 做补位；只有普通调度中进入过的组更新。

FillIdleEntries:
  使用 remaining，按 group index 从小到大访问此前未进入的组；
  普通进入组和补位进入组的 EarlyStop 行更新，其他组不更新。
```

例如普通调度第一轮进入 4 个组、第二轮进入 2 个组，之后已经没有普通候选：

```text
used_group_entries = 6
remaining = 2
```

`FillIdleEntries` 再选两个此前未进入的最小 index 组，使总访问次数达到 8。这两个组的 EarlyStop 行允许更新。

### 5.4 `K = 0`

没有普通候选：

```text
AllGroups:
  16 组全部更新，不消耗普通 HISO/SISO entry。

EnteredGroupsOnly:
  没有任何组进入，16 组全部不更新。

FillIdleEntries:
  按 index 访问 Group 0-7；Group 0-7 更新，Group 8-15 不更新。
```

## 6. 调度器内部状态

调度函数使用一个局部数组记录本次调用中哪些组已经进入：

```cpp
std::array<bool, kLevel56GroupedGroupCount> group_entered{};
```

普通调度选中组时置为 `true`；`FillIdleEntries` 选中补位组时也置为 `true`。

不能只依赖 `planned_hiso || planned_siso` 反推组是否进入，因为补位组不会规划 HISO/SISO code。也不能依赖 `Level56ScheduleSample`，因为 sample 只在 observability 开启时存在，算法结果不能随观测开关变化。

这个数组只是调度函数内部的临时算法状态，不新增公共调度表。

## 7. Observability

当 `LEVEL56_SCHEDULE_OBSERVABILITY_ENABLE=true` 时，补位 entry 也计入：

```text
Level56ScheduleSample::group_entry_counts
Level56ScheduleSample::total_group_entries
Level56ScheduleSample::rounds
```

补位轮的 `selected_groups` 记录按 index 选择的组。由于补位不运行 HISO/SISO，相关 code 的：

```text
planned_hiso = false
planned_siso = false
assigned_entry_slot = -1
assigned_core = -1
```

算法使用局部 `group_entered`，因此 observability 开关关闭时行为完全相同。

## 8. 执行和写回

不允许更新的 EarlyStop 行被设置为 `Unscheduled`。现有执行路径会跳过 `Unscheduled`，其 `produced_rows[row]` 保持 `false`。

现有 `writeback_tile()` 仅写回 `produced=true` 的行，因此无需修改：

```text
不执行 EarlyStop action
不乘 ALPHA
不重新量化
不覆盖 tile_out
不更新 Level 6 history
```

EarlyStop action 1-8、Chase、HISO/SISO core 和 Level 1-4 均不修改。

## 9. `ofec_sweep3` 切换方式

只修改以下常量即可选择一次 sweep 使用的模式：

```cpp
static constexpr newcode::Level56EarlyStopGroupUpdateMode
    kLevel56EarlyStopGroupUpdateMode =
        newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries;
```

可替换为：

```cpp
newcode::Level56EarlyStopGroupUpdateMode::AllGroups
newcode::Level56EarlyStopGroupUpdateMode::EnteredGroupsOnly
newcode::Level56EarlyStopGroupUpdateMode::FillIdleEntries
```

`build_config()` 将该值写入 `Params`。CSV 使用字段：

```text
level56_early_stop_group_update_mode
```

对应值为：

```text
all_groups
entered_groups_only
fill_idle_entries
```

## 10. 修改范围

本功能只修改：

```text
include/newcode/params.hpp
src/rx/ofec/detail/ofec_level56_shared.ipp
apps/ofec_sweep3.cpp
apps/ofec_level56_shared_regression_check.cpp
doc/ofec_level56_unentered_group_early_stop_no_update_plan.md
```

不修改公共 sweep scenario 生成逻辑，不让一个运行自动展开三种模式。

## 11. 回归验收

回归测试覆盖：

1. `AllGroups` 保持历史行为。
2. `EnteredGroupsOnly` 在 `K=0` 和 `K>8` 时抑制未进入组。
3. `FillIdleEntries` 在 `K=0` 时只访问 Group 0-7。
4. 普通调度使用 6 次 entry 后，补位选择两个最小 index 的未进入组。
5. 普通调度已经用满 8 次时，`FillIdleEntries` 与 `EnteredGroupsOnly` 一致。
6. observability 关闭不改变调度动作结果。
7. 被抑制行不产生输出，也不调用 SISO core。
