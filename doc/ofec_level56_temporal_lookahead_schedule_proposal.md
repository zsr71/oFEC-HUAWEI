# oFEC 第五/六级共享时间窗口调度新方案

## 1. 方案概述

当前第五级和第六级共享模式下，每个解码时刻处理一批共 64 个 code。新方案计划把调度视图扩展到相邻的三个解码时刻：在当前解码时刻 `t=1`，同时接收 `t=0`、`t=1` 和 `t=2` 三批 code。

```text
t=0：64 个 code，上一时刻已经处理过，可能包含尚未完成解码的遗留 code
t=1：64 个 code，当前时刻需要解码
t=2：64 个 code，未来时刻即将解码

调度输入总量：64 + 64 + 64 = 192 个 code
```

该方案的目标是利用已经解码的历史 code、当前待解码 code 以及未来 code 的早停信息，扩大调度器的观察范围，为后续的资源分配和优先级决策提供更多信息。

## 2. 调度信息划分

192 个 code 在调度阶段按顺序划分为两类信息：

```text
前 64 个 code（对应 t=0）：获取完整的解码情况
后 128 个 code（对应 t=1、t=2）：获取早停情况
```

其中：

- `t=0` 的 64 个 code 已经在前一时刻处理过，调度器可以获得其完整解码状态，其中可能包含应当解码但尚未解码的遗留 code；
- `t=1` 的 64 个 code 是当前需要处理的 code，调度器先获得早停判断信息；
- `t=2` 的 64 个 code 是未来即将处理的 code，调度器先获得早停判断信息，不要求此时完成完整解码。

这里的“完整解码情况”和“早停情况”是调度所使用的信息口径，不等同于 192 个 code 在当前时刻全部执行完整解码。

## 3. 调度 entry 和索引映射

### 3.1 调度 entry

第一版中，每个 code 使用一个统一的调度 entry 表示：

```cpp
struct Level56TemporalEntry {
  std::size_t shared_index;     // 192-code 调度视图中的统一索引：0..191
  int time_offset;              // -1、0、+1，分别表示 t=0、t=1、t=2
  std::size_t source_level;     // 来源级别：5 或 6
  std::size_t local_code_index; // 来源级别内的 code 索引：0..31

  InfoType info_type;           // DecodeInfo 或 EarlyStopInfo
  CompletedDecodeInfo decode_info; // t=0 的完整解码情况
  bool early_stop_hit;          // 逐行早停结果；随 code 状态滚动到 t=0
};

struct CompletedDecodeInfo {
  DecodeStatus status; // 解码状态
  DecodeAction action; // 实际采用的解码动作
  bool produced;       // 是否真正产生了解码输出
};
```

其中：

- `t=0` 的 entry 使用 `decode_info` 保存已完成批次的解码情况；
- `t=1`、`t=2` 通过逐行判定生成 `early_stop_hit`；该标志随后随 code 一起滚动到 `t=0`，用于区分 EarlyStop code 和真正需要 HISO/SISO 的遗留 code；
- `info_type` 用来明确当前 entry 携带的是完整解码信息还是早停信息；
- `status` 表示该 code 的解码状态；
- `action` 表示该 code 实际采用的解码动作，例如 early-stop、HISO、SISO 或未调度；
- `produced` 表示该 code 是否真正产生了解码输出。

第一版中，`t=0` 不额外保存解码质量、资源占用、完成时间、`output_valid` 或 `writeback_done`。完整的 LLR、`lout` 和历史信息也不复制到调度 entry 中，继续由现有解码数据结构保存。

### 3.2 逐行早停判定

`t=1` 和 `t=2` 的早停信息先采用与当前第五/六级共享方案一致的逐行判定方式。对每一个 code 独立执行 early-stop 条件判断，并直接得到一个布尔状态：

```text
early_stop_hit = true   -> 该 code 满足早停条件
early_stop_hit = false  -> 该 code 不满足早停条件
```

本方案固定采用 `Level56EarlyStopGroupUpdateMode::AllGroups` 的 EarlyStop 执行语义。对于当前参与处理的 `t=1` 批次，只要某个 code 满足：

```text
early_stop_hit = true
```

该 code 就执行 `EarlyStopAction`，不要求其所在四行组获得普通 group entry，也不需要额外分配 EarlyStop 补位 entry。`t=2` 在本次窗口中只提供逐行 EarlyStop 判断结果，不执行 action；当它滚动为下一窗口的 `t=1` 后，再按照上述 `AllGroups` 规则执行。

第一版暂不引入 `bind group`，不做跨行或跨级的早停状态绑定，也不区分 raw/effective 两层结果。`t=1` 的逐行结果用于当前调度；`t=2` 的逐行结果作为未来调度信息。由于四行缓冲区隔离了相邻时刻的输入，`t=2` 滚动为下一窗口的 `t=1` 后沿用该逐行 EarlyStop 结果。

上面的类型名是方案阶段的语义表示，正式落地时再结合现有代码类型确定实际定义。

### 3.3 固定线性布局

192 个 code 采用固定线性布局：先按 `t=0`、`t=1`、`t=2` 排列，每个时刻内部固定先放第五级的 32 个 code，再放第六级的 32 个 code。

```text
shared[0..31]       = t=0，Level 5 code[0..31]
shared[32..63]      = t=0，Level 6 code[0..31]

shared[64..95]      = t=1，Level 5 code[0..31]
shared[96..127]     = t=1，Level 6 code[0..31]

shared[128..159]    = t=2，Level 5 code[0..31]
shared[160..191]    = t=2，Level 6 code[0..31]
```

设 `time_index` 分别取 `0`、`1`、`2`，对应 `t=0`、`t=1`、`t=2`，则统一索引的计算方式为：

```text
time_base = time_index * 64

Level 5：shared_index = time_base + local_code_index
Level 6：shared_index = time_base + 32 + local_code_index

其中 local_code_index 的范围为 0..31。
```

例如：

```text
shared[0]   = t=0，Level 5，local code 0
shared[32]  = t=0，Level 6，local code 0
shared[64]  = t=1，Level 5，local code 0
shared[96]  = t=1，Level 6，local code 0
shared[128] = t=2，Level 5，local code 0
shared[160] = t=2，Level 6，local code 0
```

该布局负责建立稳定的输入、状态和结果映射。它本身不表示调度优先级，也不预先规定资源分配或跨时刻抢占规则。

## 4. 当前确定的调度规则

调度时先检查 `t=0` 的 64 个 code 中，是否存在“应当解码但是没有解码”的遗留 code。

本方案中的 `X`、`K1`、`K2` 沿用原第五/六级四行分组共享调度中 `K` 的统计口径。每个时刻的 64 个 code 固定分为 16 个四行组：

```text
K  = 第一轮调度开始前，包含普通非 EarlyStop 候选的非空组数量
X  = t=0 中仍包含应解码但未解码 code 的遗留组数量
K1 = t=1 中包含普通非 EarlyStop 候选的非空组数量
K2 = t=2 中包含普通非 EarlyStop 候选的非空组数量
```

原共享方案每个时刻最多使用 8 次 group entry，因此 `t=1` 和 `t=2` 两个时刻合计对应：

```text
8 + 8 = 16 次 group entry
```

后面的 `16 - X` 只用于组级负载判断，不表示精确的剩余 group-entry 数量。`X/K1/K2` 都是各自时刻第一轮开始前的非空组数；同一个组在原多轮调度中仍可能重复进入。

### 4.1 `t=0` 没有遗留 code

如果 `t=0` 不存在应当解码但没有解码的 code，则直接处理 `t=1`。其中普通非 EarlyStop code 仍复用原第五/六级 Group4 调度规则，EarlyStop code 则固定采用 `AllGroups` 语义。

```text
t=0 无遗留 code
    -> 按原 Group4 规则调度 t=1 的普通非 EarlyStop code
    -> t=1 中所有 EarlyStop 命中 code 执行 EarlyStopAction
```

### 4.2 `t=0` 存在遗留 code

如果 `t=0` 存在应当解码但没有解码的 code，并且这些 code 对应 `X` 个遗留组，则进一步读取 `t=1` 和 `t=2` 对应的 `K1`、`K2`：

```text
如果 K1 + K2 < 16 - X：
    先补解 t=0 中应当解码但没有解码的 code
    t=0 补解完成后，剩余轮次继续解码 t=1

如果 K1 + K2 >= 16 - X：
    优先解码 t=1 的 code
```

因此，等号情况归入第二个分支，即优先解码 `t=1`。`t=2` 的信息通过 `K2` 参与本次决策，但不会在本次调度中直接成为解码对象。

在该分支中，`t=0` 的遗留 code 保持不动：本次不为它们分配 HISO/SISO entry，不执行新的解码动作，也不覆盖原有状态。它们继续保持原来的：

```text
final_action = Unscheduled
produced = false
```

这些遗留 code 在本次窗口中保持不动，随后随旧 `t=0` 一起移出三时刻窗口，不进入下一窗口，也不再继续参与后续调度。

第一个判断条件也可以等价写成：

```text
X + K1 + K2 < 16
```

它表示把 `t=0` 的遗留组和 `t=1/t=2` 的普通非 EarlyStop 非空组放在一起观察时，总组级负载低于预设阈值 16。该阈值只用于决定优先补解 `t=0` 还是优先处理 `t=1`，不用于精确计算实际需要的 group-entry 数量。

### 4.3 补解 `t=0` 时的两阶段调度

当判断结果为先补解 `t=0` 时，本次解码时刻的 8 次 group-entry 预算按照两个阶段依次使用。

第一阶段处理 `t=0`：

1. 调度候选只包含 `t=0` 中应当解码但未解码的 code。
2. 保留原来的固定四行分组，不重新改变 code 与组的对应关系。
3. 组进入顺序复用原 Group4 load-sorted multiround 规则：按照当前剩余未调度 code 数量从高到低选择，负载相同时按原组号顺序处理。
4. 组内 code 选择复用原规则：SISO 优先选择 `SisoOnly` 候选，否则选择 `HisoOrSiso` 候选；HISO 再从剩余的 `HisoOrSiso` 候选中选择一个不同的 code。
5. 同一遗留组仍有未调度 code 时，可以按照原多轮规则再次进入；已经获得 entry slot 的 code 不得重复选择。

第二阶段处理 `t=1`：

```text
t=0 已经没有遗留候选，并且本次 8 次 group-entry 预算尚未用完
    -> 将剩余 group-entry slot 交给 t=1
    -> t=1 继续按照原 Group4 load-sorted multiround 规则调度
```

这里的“剩余轮次”具体表现为本次尚未使用的 group-entry slot。设第一阶段处理 `t=0` 实际使用了 `used_t0` 次 group entry，则 `t=1` 最多继续使用：

```text
remaining_for_t1 = 8 - used_t0
```

如果 `t=0` 的补解已经用完 8 次 group entry，则本次不再为 `t=1` 安排普通 HISO/SISO 解码；`t=1` 中命中 EarlyStop 的 code 仍执行 `EarlyStopAction`。

如果 `t=0` 的遗留候选提前耗尽，则只将尚未使用的 entry 交给 `t=1`；如果 8 次 entry 已经用完，则本次不再为 `t=1` 的普通非 EarlyStop code 分配资源。

#### 4.3.1 `t=1` 独立重新开始 Group4 调度

`t=0` 补解结束后，`t=1` 不延续 `t=0` 的组排序、剩余负载表或多轮状态，而是根据 `t=1` 自己的 64 个 code 独立建立一套 Group4 调度状态：

```text
t=1 的 64 个 code
    -> 固定划分为 16 个四行组
    -> 排除已经命中 EarlyStop 的 code
    -> 对普通非 EarlyStop code 执行 Hybrid classify-only
    -> 统计每组的初始待调度 code 数量 initial_count_i
    -> 统计 K1 = count(initial_count_i > 0)
```

`t=0` 和 `t=1` 的组号都可以是 `0..15`，但属于两套彼此独立的时间批次。`t=0` 的 Group 0 不会和 `t=1` 的 Group 0 合并，也不会共同参与同一次组负载排序。

#### 4.3.2 `t=1` 的可用预算

设 `t=0` 补解实际使用：

```text
used_t0 次 group entry
```

则 `t=1` 本次可使用的总预算为：

```text
B = 8 - used_t0
```

`B` 同时约束 `t=1` 的第一轮和所有后续轮次：

```text
t=1_total_group_entries <= B
```

资源槽连续使用：`t=0` 使用 entry slot `0..used_t0-1`，`t=1` 使用剩余的 entry slot `used_t0..7`。如果 `B == 0`，`t=1` 本次不安排普通 HISO/SISO 解码，但其中所有 `early_stop_hit=true` 的 code 仍按照 `AllGroups` 规则执行 `EarlyStopAction`。EarlyStopAction 不占用这里的 HISO/SISO group-entry 预算。

#### 4.3.3 按 `K1` 和 `B` 执行预算化多轮调度

`t=1` 的 Group4 调度保持原方案的分支逻辑，只把原来的最大预算 8 替换为本次剩余预算 `B`。

这是 `t=1` 调度与原第五/六级 Group4 load-sorted multiround 方案之间唯一的规则差异：

```text
原方案：max_group_entries = 8
本分支：max_group_entries = B = 8 - used_t0
```

除 `max_group_entries` 外，普通 HISO/SISO 调度的以下行为全部与原方案保持一致，不增加新的排序或选择规则：

```text
64 个 code 固定分成 16 个四行组
initial_count_i 和 remaining_unscheduled_count_i 的定义
K1 的非空组统计口径
按组负载从高到低排序
组负载相同时按原组号从小到大排序
SisoOnly 和 HisoOrSiso 的资源资格
组内 SISO/HISO code 的选择顺序
同一组的多轮重入规则
已预约 code 不得重复选择
未调度 code 的 Unscheduled / produced=false 语义
```

EarlyStopAction 独立固定采用 `AllGroups` 语义，不受 `t=1` 的预算 `B`、普通调度进入组集合或多轮调度结果影响。

实现时应优先复用原 Group4 调度函数，并把最大 group-entry 数量参数化为 `max_group_entries`；不应为 `t=1` 复制出另一套行为不同的调度逻辑。

当 `K1 > B` 时：

```text
按 initial_count_i 从高到低排序
负载相同时按原组号从小到大排序
只选择前 B 个非空组
每个选中组进入一次
预算用完，不再进行后续轮次
```

当 `K1 == B` 时：

```text
全部 K1 个非空组各进入一次
预算恰好用完
不再进行后续轮次
```

当 `0 < K1 < B` 时：

```text
第一轮让全部 K1 个非空组各进入一次
重新统计每组 remaining_unscheduled_count_i
如果仍有预算和未调度候选，则进入下一轮
后续每轮按 remaining_unscheduled_count_i 从高到低重新排序
直到 t=1 累计使用 B 次 entry，或所有待调度候选耗尽
```

当 `K1 == 0` 时，`t=1` 没有普通非 EarlyStop 候选，不使用 HISO/SISO entry。

每个被选中的 `t=1` 组继续复用原方案的组内规则：

```text
SISO：优先选择一个 SisoOnly；没有时选择一个 HisoOrSiso
HISO：从剩余 HisoOrSiso 中选择最多一个不同的 code
已经获得 entry slot 的 code 不得在后续轮次重复选择
```

全部 `t=1` 调度轮次结束后，仍未获得 HISO/SISO entry 的普通非 EarlyStop code 保持 `Unscheduled / produced=false`，并在下一窗口顺延为新的 `t=0` 遗留候选。

### 4.4 `t=1` 的顺延

无论是因为本次 entry 预算已经用完，还是因为 `t=1` 只完成了部分调度，所有未在本次获得 HISO/SISO entry slot 的 `t=1` 非 EarlyStop code 都保持：

```text
final_action = Unscheduled
produced = false
```

到下一个解码时刻，原来的 `t=1` 顺延为新的 `t=0`。其中仍应解码但未解码的 code 成为新的历史遗留候选，并按照相同方式统计新的 `X`。

这里的顺延只适用于当前 `t=1`，不适用于已经从窗口中移出的旧 `t=0`。旧 `t=0` 未补解的遗留 code 在窗口移出后不再保留。

`t=1` 中命中 EarlyStop 的 code 不参与普通 HISO/SISO 候选选择，并统一按照 `AllGroups` 规则执行 `EarlyStopAction`：

```text
early_stop_hit = true
    -> 执行 EarlyStopAction
    -> 不检查所在组是否获得普通 group entry
    -> 不使用空闲 entry 为 EarlyStop 单独补位
```

因此，即使所在组从未获得普通调度机会，EarlyStop code 仍然执行 action。EarlyStop code 本身不占用 HISO/SISO core，也不减少 `t=0` 或 `t=1` 可使用的普通 group-entry 数量。

### 4.5 时间窗口滚动与 `X/K1/K2` 更新

每次调度完成后，三个时刻的窗口整体向前滚动一个解码时刻：

```text
本次窗口：       [旧 t=0] [旧 t=1] [旧 t=2]
本次调度：          历史      当前      未来
                     |         |         |
下一次窗口：       [移出] [新 t=0] [新 t=1] [新 t=2]
                              |         |         |
                           旧 t=1     旧 t=2     新 t=3
```

具体对应关系为：

```text
旧 t=0 -> 从三时刻调度窗口中移出
旧 t=1 -> 成为下一次窗口的新 t=0
旧 t=2 -> 成为下一次窗口的新 t=1
新 t=3 -> 进入下一次窗口，成为新 t=2
```

下一次调度所使用的 `X/K1/K2` 按滚动后的三个时刻重新统计：

```text
X_new  <- 旧 t=1 中遗留的应解码但未解码组
K1_new <- 旧 t=2 中包含普通非 EarlyStop 候选的非空组
K2_new <- 新 t=3 中包含普通非 EarlyStop 候选的非空组
```

旧 `t=1` 中已经完成解码的 code，滚动为新 `t=0` 后保留其 `DecodeStatus / DecodeAction / produced` 和 `early_stop_hit`；旧 `t=1` 中仍未获得 HISO/SISO 解码机会的普通 code，滚动后继续保持 `Unscheduled / produced=false`。其中只有 `early_stop_hit=false` 的普通 code 参与 `X_new` 的统计；`early_stop_hit=true` 的 code 已经按照 `AllGroups` 规则执行 `EarlyStopAction`，不计入 `X_new`。

旧 `t=2` 滚动为新 `t=1` 后，沿用其逐行 EarlyStop 判断结果并统计 `K1_new`。由于四行缓冲区已经隔离相邻时刻的输入，旧 `t=1` 的解码写回不会改变旧 `t=2` 的输入。新进入的 `t=3` 作为新 `t=2`，执行逐行 EarlyStop 判断并统计 `K2_new`。

### 4.6 遗留 code 的 Hybrid 分类

`t=0` 遗留 code 再次参与 Group4 调度时，重新执行一次 Hybrid classify-only，并根据本次输入生成新的 `hybrid_class` 和资源资格。

在当前方案的理想条件下，未调度 code 没有产生新的写回，且分类参数保持不变，因此重新分类的结果应与上一次分类结果一致。采用重新分类主要是为了复用现有的输入准备和分类流程，避免在时间窗口滚动时额外保存和维护遗留 code 的 `hybrid_class`、资源资格等调度上下文。具体的重算时机属于实现阶段问题。

重新分类不改变该 code 原有的解码状态：如果它此前未获得解码机会，仍然保持：

```text
final_action = Unscheduled
produced = false
```

重新分类只更新它重新参加 Group4 调度所需的候选属性。

### 4.7 遗留 code 的判定

逐 code 定义：

```text
is_pending_decode =
    early_stop_hit == false
    && final_action == Unscheduled
    && produced == false
```

只有满足 `is_pending_decode == true` 的 code 才属于应当解码但未解码的遗留 code。一个四行组中只要至少有一个 code 满足该条件，该组就计入 `X`。

以下 code 不属于“应当进入 HISO/SISO 解码但未解码”的遗留 code：

```text
已经按照 AllGroups 规则执行 EarlyStopAction 的 EarlyStop code
已经 HISO/SISO 解码并产生输出的 code
```

整体判断可以写成：

```text
检查 t=0 是否存在应解码但未解码的 code
    |
    +-- 否 -> 按原 Group4 规则调度 t=1 的普通 code
    |        t=1 的 EarlyStop code 按 AllGroups 执行 action
    |
    +-- 是，对应 X 组
            |
            +-- K1 + K2 < 16 - X  -> 先补解 t=0，剩余 entry 解码 t=1
            |
            +-- K1 + K2 >= 16 - X -> 优先解码 t=1
```

## 5. 与当前五、六级共享流程的关系

新方案保留当前五、六级共享的基本边界：每个解码时刻仍由第五级和第六级共同形成 64 个 code 的共享批次；变化只在于调度器一次观察三个连续时刻的批次。

可以将当前流程理解为：

```text
当前方案：一次观察 1 个 64-code 批次
新方案：一次观察 3 个连续的 64-code 批次，共 192 个 code
```

每个时刻内部的普通 HISO/SISO 调度仍沿用原方案的固定四行分组、非空组统计和 group-entry 概念。EarlyStop 更新则固定采用 `AllGroups`：当前 `t=1` 中所有命中 EarlyStop 的 code 均执行 `EarlyStopAction`，不受所在组是否进入及普通 entry 预算影响。新方案在其外层增加跨时刻判断，用 `X`、`K1`、`K2` 决定本次优先补解 `t=0`，还是优先处理 `t=1`。输入准备、第五/六级来源映射以及最终结果写回仍需要能够区分所属的解码时刻。

## 6. 边界情况

序列开始和结束时只会缺少一侧的批次，不考虑历史批次和未来批次同时缺失的情况。

### 6.1 开始阶段缺少历史批次

开始阶段不存在历史 `t=0` 批次时：

```text
X = 0
```

不生成虚假的历史遗留 code，也不为缺失的历史批次分配资源。第一个实际批次直接作为当前 `t=1`，如果未来批次存在则同时作为 `t=2` 提供 lookahead 信息。当前批次的普通非 EarlyStop code 按原 Group4 规则调度，EarlyStop code 固定采用 `AllGroups` 语义：

```text
t=0 无遗留 code
    -> 按原 Group4 规则调度当前 t=1 的普通非 EarlyStop code
    -> 当前 t=1 中所有 EarlyStop 命中 code 执行 EarlyStopAction
```

### 6.2 结束阶段缺少未来批次

结束阶段不存在未来 `t=2` 批次时，不构造虚假的未来 code、EarlyStop 状态或 `K2`。直接处理当前 `t=1`：

```text
缺少未来 t=2
    -> 只处理当前 t=1
    -> 不执行 t=2 的 lookahead 调度
```

该情况不增加额外的 flush 或补解流程。最后一块数据在 BER 统计时会被裁掉，因此不会因为缺少未来批次而改变最终 BER 统计口径。

## 7. 后续需要讨论和确认的内容

当前没有尚未确认的方案级流程或调度规则。

EarlyStop 的组更新模式已经确定为：

```text
Level56EarlyStopGroupUpdateMode::AllGroups
```

当前 `t=1` 中所有 `early_stop_hit=true` 的 code 均执行 `EarlyStopAction`，不依赖所在组是否获得普通 group entry，也不使用空闲 entry 执行补位。该规则与普通 HISO/SISO 的预算分配相互独立。

此前关于 `t=0` 遗留 code 状态口径的问题已经确定：`CompletedDecodeInfo` 仍只保存 `DecodeStatus / DecodeAction / produced` 三个字段；顶层 entry 的 `early_stop_hit` 独立随 code 状态滚动到 `t=0`，不计入这三个完整解码结果字段。遗留 code 继续使用以下条件判定：

```text
is_pending_decode =
    early_stop_hit == false
    && final_action == Unscheduled
    && produced == false
```

缓存组织、具体状态字段搬运、吞吐率统计和回归验证方式属于实现阶段内容，待实际修改代码时再单独确定。
