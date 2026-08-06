# oFEC 第五/六级共享时间窗口调度新方案

## 1. 方案概述

当前第五级和第六级共享模式下，每个解码时刻处理一批共 64 个 code。新方案计划把调度视图扩展到相邻的三个解码时刻：在当前解码时刻 `t=1`，同时接收 `t=0`、`t=1` 和 `t=2` 三批 code。

```text
t=0：64 个 code，已经完成上一时刻解码
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

- `t=0` 的 64 个 code 已经在前一时刻完成解码，调度器可以获得其完整解码结果及相关状态；
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
  bool early_stop_hit;          // t=1、t=2 的逐行早停结果
};

struct CompletedDecodeInfo {
  DecodeStatus status; // 解码状态
  DecodeAction action; // 实际采用的解码动作
  bool produced;       // 是否真正产生了解码输出
};
```

其中：

- `t=0` 的 entry 使用 `decode_info` 保存已完成批次的解码情况；
- `t=1`、`t=2` 的 entry 使用 `early_stop_hit`，保存对应 code 的逐行早停结果；
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

第一版暂不引入 `bind group`，不做跨行或跨级的早停状态绑定，也不区分 raw/effective 两层结果。`t=1` 的逐行结果用于当前调度；`t=2` 的逐行结果只作为未来调度的预先信息，真正处理 `t=2` 时再重新执行该 code 的 early-stop 判定。

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

后面的 `16 - X` 表示在考虑 `t=0` 的 `X` 个遗留组后，两个时刻合计 16 次 group entry 中剩余的组级容量。

### 4.1 `t=0` 没有遗留 code

如果 `t=0` 不存在应当解码但没有解码的 code，则保持原第五/六级共享方案的行为，直接解码 `t=1` 的部分。

```text
t=0 无遗留 code
    -> 按原第五/六级共享方案调度并解码 t=1
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

第一个判断条件也可以等价写成：

```text
X + K1 + K2 < 16
```

它表示把 `t=0` 的遗留组和 `t=1/t=2` 的普通非 EarlyStop 非空组放在一起观察时，总组数小于两个时刻合计的 16 次 group-entry 容量。

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

如果 `t=0` 的补解已经用完 8 次 group entry，则本次不再解码 `t=1`。

### 4.4 `t=1` 的顺延

无论是因为本次 entry 预算已经用完，还是因为 `t=1` 只完成了部分调度，所有未在本次获得 HISO/SISO entry slot 的 `t=1` 非 EarlyStop code 都保持：

```text
final_action = Unscheduled
produced = false
```

到下一个解码时刻，原来的 `t=1` 顺延为新的 `t=0`。其中仍应解码但未解码的 code 成为新的历史遗留候选，并按照相同方式统计新的 `X`。

整体判断可以写成：

```text
检查 t=0 是否存在应解码但未解码的 code
    |
    +-- 否 -> 按原共享方案解码 t=1
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

每个时刻内部仍沿用原方案的固定四行分组、非空组统计和 group-entry 概念。新方案在其外层增加跨时刻判断，用 `X`、`K1`、`K2` 决定本次优先补解 `t=0`，还是优先处理 `t=1`。输入准备、第五/六级来源映射以及最终结果写回仍需要能够区分所属的解码时刻。

## 6. 边界情况

序列开始和结束时可能无法同时取得完整的三个时刻，例如开始时不存在已经完成的 `t=0` 批次，结束时不存在未来的 `t=2` 批次。

边界情况暂不影响主体方案的定义，后续统一处理。届时再确定采用补齐、缩短调度窗口或其他策略。

## 7. 后续需要讨论和确认的内容

目前仍需要继续讨论和确认的内容为：

1. 当判断结果为优先解码 `t=1` 时，本次未补解的 `t=0` 遗留 code 后续是否继续保留，以及允许保留多少个时刻。
2. 每次调度结束后，三个时刻的 entry 和 `DecodeStatus / DecodeAction / produced` 状态如何具体滚动、更新和重新计算 `X/K1/K2`。
3. `t=0` 遗留 code 再次参与原 Group4 调度时，所需的 Hybrid 分类和资源资格是保留上一时刻结果，还是基于当前输入重新计算。
4. 最后统一处理解码序列开始和结束时缺少历史批次或未来批次的边界情况。
5. 方案落地时的缓存规模、流水线时序、吞吐率以及统计验证方式。
