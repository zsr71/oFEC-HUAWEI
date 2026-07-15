# oFEC 第五/六级共享 HISO/SISO 解码流程设计

## 1. 文档目标

本文定义一种单流 oFEC 第五级和第六级合并解码方案：

- 第五级输入 32 个 code；
- 第六级输入 32 个 code；
- 两级合并为一个 64-code shared batch；
- 64 个 code 共享一组 HISO 核和一组 SISO 核；
- early-stop、hybrid 分类和跨级优先级在 64-code 域内统一处理；
- HISO/SISO 解码完成后，再按来源级别分别进行后处理和写回。

本文只定义解码流程和模块职责，不修改当前代码。

本文使用 `HISO` 表示共享硬输入、软输出（或当前方案中承担 hard-finish 的）解码资源。当前代码中的资源参数仍命名为 `HIHO_ACTIVE_LIST`，正式实现前需要再统一术语，避免设计名称和代码名称混淆。

## 2. 设计原则

第五/六级共享模式遵循以下原则：

1. 合并时不交错排序。共享数组固定先放第五级 32 个 code，再放第六级 32 个 code。
2. early-stop 对 64 个 code 分别判断，但每个 code 使用所属级别的参数。
3. 当 `EARLY_STOP_BIND_GROUP_SIZE=4` 时，group binding 只允许发生在同一级内部。
4. hybrid 阶段只分类，不翻 bit、不执行 BCH hard decode、不生成最终 hard 输出。
5. 跨第五/六级的资源竞争、优先级和最终路径选择全部放在统一优先级处理阶段。
6. HISO MUX 和 SISO MUX 只负责物理或逻辑路由，不负责优先级仲裁。
7. HISO/SISO 解码后，第五级和第六级分别使用各自的 alpha、beta、early-stop action 和地址映射完成写回。

## 3. 64-code 合并顺序

第五级和第六级分别完成输入准备后，组成一个 64-code shared batch：

```text
shared[0..31]  = Level 5 code[0..31]
shared[32..63] = Level 6 code[0..31]
```

第一版不在合并阶段使用以下交错顺序：

```text
Level5 row0, Level6 row0, Level5 row1, Level6 row1, ...
```

原因是合并顺序只用于保存来源和建立稳定索引，不应该隐式表达资源优先级。交错、公平、第五级优先和第六级优先都属于后面的统一优先级策略。

每个 shared entry 至少保存：

```cpp
struct Level56SharedCodeEntry {
  std::size_t shared_row;       // 0..63
  std::size_t source_level;     // 5 或 6
  std::size_t source_local_row; // 0..31
  std::size_t source_global_row;

  LinVector lin;
  LchVector lch;

  bool early_stop_hit;
  HybridRowClass hybrid_class;
  ResourceEligibility eligibility;
  FinalAction final_action;
};
```

`source_level` 和 `source_local_row` 在整个流程中不能丢失，因为后续需要据此选择对应参数并恢复第五/六级写回映射。

## 4. 总体流程

```text
前四级按现有流程处理并写回
                |
                v
   准备 Level 5 的 32 个 code
   准备 Level 6 的 32 个 code
                |
                v
     合并为 64-code shared batch
   [L5 0..31][L6 0..31]，不交错
                |
                v
        64 个 code 做 early-stop
                |
                v
   未早停 code 做 hybrid classify-only
                |
                v
       建立 64 行分类和资格状态表
                |
                v
       跨第五/六级统一优先级处理
   - 决定 HISO / SISO / Unscheduled
   - 支持公平 / L5 优先 / L6 优先
                |
          +-----+-----+
          |           |
          v           v
      HISO MUX     SISO MUX
      只做路由      只做路由
          |           |
          v           v
    共享 HISO 核   共享 SISO 核
          |           |
          +-----+-----+
                |
                v
      合并 64 行最终执行结果
                |
          +-----+-----+
          |           |
          v           v
    拆回 Level 5  拆回 Level 6
    独立后处理      独立后处理
    独立写回        独立写回
```

## 5. 第一步：分别准备第五/六级输入

第五级和第六级必须分别执行 `prepare_tile_inputs()`，不能先把原始 tile 矩阵拼接后只准备一次。

原因包括：

- 两级的 `tile_top_row_global` 不同；
- 两级的历史信息读取地址不同；
- 两级可能使用不同的 alpha、beta 和 early-stop 参数；
- 两级的 trace、row lookup 和写回映射不同。

准备阶段产生：

```text
Level5Prepared: 32 x 256 Lin/Lch + Level5 row lookup
Level6Prepared: 32 x 256 Lin/Lch + Level6 row lookup
```

然后只在 shared 调度视图中把两组 prepared rows 合并为 64 行。原始的两个 `TilePrepared` 继续保留，供最终分级执行和写回使用。

第一版共享模式建议限制：

```text
TILES_PER_WIN == 6
CHASE_SBR == 2
TILE_OVERLAP_BR == 0
Level 5 和 Level 6 都是 soft/hybrid tile
```

这些限制可以避免第六级输入依赖第五级即时写回、hard-tile history input 或 overlap 语义，从而先验证共享资源本身对 BER 的影响。

## 6. 第二步：64 个 code 做 early-stop

64 行分别执行现有 early-stop 条件：

```text
Level 5 code 使用 Level 5 的 early-stop 参数
Level 6 code 使用 Level 6 的 early-stop 参数
```

输出分成：

```text
EarlyStopAction
NotEarlyStop
```

early-stop 命中的 code：

- 不进入 hybrid 分类；
- 不参与跨级优先级计算；
- 不进入 HISO MUX；
- 不进入 SISO MUX；
- 不占用 HISO/SISO 核；
- 在执行结果合并阶段，使用所属级别的 early-stop action 生成输出。

### 6.1 `EARLY_STOP_BIND_GROUP_SIZE=4`

绑定严格限制在同一级内部：

```text
Level 5 groups:
  L5[0..3], L5[4..7], ..., L5[28..31]

Level 6 groups:
  L6[0..3], L6[4..7], ..., L6[28..31]
```

禁止因为 shared batch 的存在而跨级绑定，例如：

```text
L5[30], L5[31], L6[0], L6[1]
```

也禁止按照公平交错次序形成跨级 group。group binding 是 early-stop 语义，不是资源公平策略。

## 7. 第三步：hybrid 只做分类

所有未 early-stop 的 code 进入 hybrid classify-only。

分类阶段输入：

```text
Lin[256]
所属级别的 classifier 参数
```

分类阶段输出：

```text
Clean
ParityOnly
OneMain
OneMainPlusParity
TwoMain
Suspicious
HardFail
```

分类阶段禁止执行：

- 禁止翻转 hard-decision bit；
- 禁止调用 BCH HISO/HIHO decoder；
- 禁止重算 overall parity 后直接生成结果；
- 禁止生成 `y2`；
- 禁止提前占用 HISO 核；
- 禁止提前调用 Chase SISO。

这一阶段只回答“这个 code 具备走哪些路径的资格”，不回答“最后分配到哪个核”。

## 8. 第四步：建立资源资格

分类结果转换成资源资格。资格和最终动作必须使用两个不同的字段。

```cpp
enum class ResourceEligibility {
  None,
  FreeFinish,
  HisoOnly,
  SisoOnly,
  HisoOrSiso,
};

enum class FinalAction {
  EarlyStopAction,
  FreeFinish,
  HisoDecode,
  SisoDecode,
  Unscheduled,
};
```

推荐的基础映射如下：

| 输入状态或分类        | 资源资格       | 说明                             |
| --------------------- | -------------- | -------------------------------- |
| early-stop 命中       | `None`       | 最终为`EarlyStopAction`        |
| `Clean`             | `FreeFinish` | 不占 HISO/SISO，使用免费完成路径 |
| `ParityOnly`        | `HisoOrSiso` | 可被回收到 HISO，也可保留给 SISO |
| `OneMain`           | `HisoOrSiso` | 可被回收到 HISO，也可保留给 SISO |
| `OneMainPlusParity` | `HisoOrSiso` | 可被回收到 HISO，也可保留给 SISO |
| `TwoMain`           | `HisoOrSiso` | 可被回收到 HISO，也可保留给 SISO |
| `Suspicious`        | `SisoOnly`   | 只允许进入 SISO                  |
| `HardFail`          | `SisoOnly`   | 只允许进入 SISO                  |

如果后续分类器定义出严格的 `HisoOnly` 类别，可以直接扩展优先级处理，不需要改变 MUX 职责。

## 9. 第五步：跨第五/六级统一优先级处理

统一优先级处理是整个 shared 方案唯一负责资源仲裁的模块。

它负责：

1. 合并第五/六级的 HISO/SISO 候选；
2. 根据分类优先级生成有序候选；
3. 在 `shared_hiso_capacity` 和 `shared_siso_capacity` 约束下分配最终路径；
4. 决定每个 code 的 `FinalAction`；
5. 生成已经排好序的 HISO route request 和 SISO route request；
6. 把无法获得任一资源的 code 标记为 `Unscheduled`。

它不负责：

- 不执行 HISO decode；
- 不执行 SISO decode；
- 不负责 code-to-core 的拓扑连线；
- 不执行第五/六级写回。

### 9.1 分类优先级

当使用现有 `ParityOneAndTwoErrorPriority` 时，HISO reclaim 的类内顺序保持为：

```text
ParityOnly
  -> OneMain
  -> OneMainPlusParity
  -> TwoMain
```

`SisoOnly` code 必须保留在 SISO 候选集合中，不能为了填满 HISO 而错误地送入 HISO。

### 9.2 三种跨级优先级策略

当两个 code 的分类优先级相同时，再使用跨级策略决定第五级和第六级的先后顺序。

#### 策略 A：公平模式

公平模式在统一优先级处理内部按两级轮转取 code：

```text
L5 first candidate
L6 first candidate
L5 second candidate
L6 second candidate
...
```

注意：这里只是优先级处理产生的有序请求，不改变原始 shared batch 的 `[L5 0..31][L6 0..31]` 存储顺序。

如果某一级在当前优先级桶中没有候选，另一级可以使用剩余资源，不能为了形式公平而让核空闲。

为了避免每个 invocation 都固定从第五级开始，可选地按 shared invocation 轮换起始级别：

```text
偶数 invocation：L5 -> L6 -> L5 -> L6
奇数 invocation：L6 -> L5 -> L6 -> L5
```

第一版也可以固定 L5 起始，但 CSV 中必须记录起始级别，便于判断是否形成长期偏置。

#### 策略 B：第五级优先

在相同分类优先级下：

```text
先选择第五级候选，再选择第六级候选。
```

第五级没有候选或其候选全部已分配后，第六级使用剩余资源。

#### 策略 C：第六级优先

在相同分类优先级下：

```text
先选择第六级候选，再选择第五级候选。
```

第六级没有候选或其候选全部已分配后，第五级使用剩余资源。

建议配置枚举为：

```cpp
enum class Level56PriorityMode {
  Fair = 0,
  Level5First = 1,
  Level6First = 2,
};
```

### 9.3 HISO/SISO 路径分配

设：

```text
H = shared_hiso_capacity
S = shared_siso_capacity
```

统一优先级处理首先统计：

```text
siso_only_count
hiso_or_siso_count
```

所有 `HisoOrSiso` 候选初始保留在 soft pool。若 soft pool 超过 SISO 容量，则按分类优先级和跨级优先级从中回收一部分进入 HISO：

```text
soft_total = siso_only_count + hiso_or_siso_count
need_hiso_reclaim = max(soft_total - S, 0)
hiso_scheduled_count = min(need_hiso_reclaim,
                           hiso_or_siso_count,
                           H)
```

然后：

1. 选中的 `HisoOrSiso` code 标记为 `HisoDecode`；
2. 未被选中的 `HisoOrSiso` code 保留在 soft pool；
3. soft pool 再按统一优先级选择最多 S 个，标记为 `SisoDecode`；
4. 剩余 code 标记为 `Unscheduled`。

如果未来加入 `HisoOnly`，它们应先进入 HISO 有序候选，再根据 HISO 容量决定 `HisoDecode` 或 `Unscheduled`。

最终由优先级处理直接输出：

```text
ordered_hiso_requests
ordered_siso_requests
unscheduled_rows
```

并满足：

```text
ordered_hiso_requests.size() <= H
ordered_siso_requests.size() <= S
```

## 10. 第六步：HISO MUX 和 SISO MUX

### 10.1 模块职责

HISO MUX 和 SISO MUX 不负责仲裁。它们接收的是统一优先级处理已经选好的、有序的 route requests。

HISO MUX 输入：

```text
ordered_hiso_requests
当前可用 HISO core 列表
HISO 拓扑和 bypass 配置
```

HISO MUX 输出：

```text
code_to_hiso_core mapping
无法路由的 requests（如果拓扑不是全连接）
```

SISO MUX 输入：

```text
ordered_siso_requests
当前可用 SISO core 列表
SISO 拓扑和 bypass 配置
```

SISO MUX 输出：

```text
code_to_siso_core mapping
无法路由的 requests（如果拓扑不是全连接）
```

MUX 不做：

- 不计算 hybrid class priority；
- 不决定第五级或第六级优先；
- 不因为看到某一类 code 而替换已选请求；
- 不扩大 HISO/SISO 调度数量；
- 不执行解码。

第一版若使用 `MUX_GROUP_G=1` 的全局池化，可以认为所有选中的请求都能路由成功。后续加入分组拓扑或 bypass edge 后，MUX 才需要报告 `waiting/unroutable`。

### 10.2 路由失败的处理

当物理拓扑导致某个已选请求无法连接到空闲核时，建议采用以下边界：

1. MUX 返回 `unroutable_requests`；
2. 优先级处理模块根据原有有序候选做补位；
3. MUX 再尝试路由；
4. 无候选可补或达到固定尝试次数后，才标记为 `Unscheduled`。

也就是说，MUX可以报告路由可行性，但最终“换谁上、谁未调度”的决定仍属于优先级处理模块。

### 10.3 与当前代码的区别

当前单 tile 软件实现中的 `run_mux_on_soft_candidates()` 同时执行了：

```text
候选收集
预算裁剪
可选优先级排序
分组或 reconfig 路由
写回 scheduled / unscheduled 状态
```

因此当前函数名虽然是 MUX，实际职责比目标硬件模块更宽。

第五/六级共享实现应在逻辑上拆成：

```text
Level56PriorityScheduler
  -> 资源仲裁、排序、路径选择、预算裁剪

Level56HisoMux / Level56SisoMux
  -> 只根据拓扑完成 code-to-core mapping
```

第一版软件可以复用现有底层 matching/schedule helper，但不能把 class priority 或 L5/L6 priority 再放回 MUX 层。

## 11. 第七步：执行共享 HISO/SISO 解码

只有 MUX 路由成功的 code 才真正占用核并执行解码。

### 11.1 HISO 执行

`HisoDecode` code 根据 `hybrid_class` 执行对应 hard executor：

```text
ParityOnly          -> 修正 overall parity
OneMain             -> 修正一个主体 bit
OneMainPlusParity   -> 修正一个主体 bit 和 overall parity
TwoMain             -> 执行 BCH t=2 hard decode
```

分类和执行结果应严格一致。如果 classify-only 认为可以 hard-finish，但 HISO executor 失败，第一版建议直接报错，以定位分类、syndrome、地址或 parity 口径的不一致。

### 11.2 SISO 执行

`SisoDecode` code 进入 Chase SISO。资源上两级共享一组核，但每个 code 必须使用所属级别的参数：

```text
Level 5 code 使用 beta[4] 及 Level 5 Params
Level 6 code 使用 beta[5] 及 Level 6 Params
```

由于当前一次 Chase batch 只携带一套 `params_for_core`，第一版软件仿真可以：

1. 先在 64-code 域内统一选择 SISO code；
2. 再将选中的 code 拆成 Level 5 和 Level 6 两个 scheduled row list；
3. 分别调用两级现有 SISO 执行函数。

只要满足下面的总资源约束，软件分两次调用仍然正确模拟共享核数量：

```text
level5_siso_scheduled + level6_siso_scheduled
    <= shared_siso_capacity
```

HISO 同理：

```text
level5_hiso_scheduled + level6_hiso_scheduled
    <= shared_hiso_capacity
```

## 12. 第八步：分级后处理和分别写回

解码完成后，根据 `source_level` 将 64 行结果拆回：

```text
Level 5 result: 32 x 256
Level 6 result: 32 x 256
```

每个 code 按 `FinalAction` 生成结果：

| `FinalAction`     | 执行结果                                      |
| ------------------- | --------------------------------------------- |
| `EarlyStopAction` | 使用所属级别的 early-stop action              |
| `FreeFinish`      | 使用 clean/free materialize                   |
| `HisoDecode`      | 使用 HISO executor 输出                       |
| `SisoDecode`      | 使用 Chase SISO 输出                          |
| `Unscheduled`     | `produced=false`，不产生新的 extrinsic 写回 |

随后两级必须分别执行：

```text
Level 5:
  normalize（若启用）
  alpha[4]
  Level 5 writeback mapping

Level 6:
  normalize（若启用）
  alpha[5]
  Level 6 writeback mapping
```

禁止把 64 行放在一起统一 normalization，也禁止使用同一个 alpha/beta 处理两级。

第一版保持现有写回先后顺序：

```text
先写回 Level 5
再写回 Level 6
```

## 13. 最终状态约束

每个 shared invocation 必须满足：

```text
early_stop_count
+ free_finish_count
+ hiso_scheduled_count
+ siso_scheduled_count
+ unscheduled_count
= 64
```

资源约束：

```text
level5_hiso_scheduled + level6_hiso_scheduled
    <= shared_hiso_capacity

level5_siso_scheduled + level6_siso_scheduled
    <= shared_siso_capacity
```

group binding 约束：

```text
所有 bind group 的 source_level 必须相同
```

来源恢复约束：

```text
每个 shared_row 必须唯一映射回
(source_level, source_local_row)
```

## 14. 建议配置项

```cpp
bool LEVEL56_SHARED_ENABLE = false;
int LEVEL56_SHARED_HISO_ACTIVE = 32;
int LEVEL56_SHARED_SISO_ACTIVE = 32;
Level56PriorityMode LEVEL56_PRIORITY_MODE = Level56PriorityMode::Fair;
bool LEVEL56_FAIR_ALTERNATE_START = true;
```

第五/六级共享模式打开后，`SISO_ACTIVE_LIST[4]`、`SISO_ACTIVE_LIST[5]`、`HIHO_ACTIVE_LIST[4]` 和 `HIHO_ACTIVE_LIST[5]` 不再分别表示两套独立资源。共享容量由专门的 `LEVEL56_SHARED_*` 参数控制，避免把“每级预算”和“跨级总预算”混为一谈。

## 15. 建议观测字段

每个 shared invocation 至少记录：

```text
window_idx
shared_invocation
priority_mode
fair_start_level

level5_early_stop
level6_early_stop
level5_free_finish
level6_free_finish

level5_hiso_candidates
level6_hiso_candidates
level5_siso_candidates
level6_siso_candidates

level5_hiso_scheduled
level6_hiso_scheduled
level5_siso_scheduled
level6_siso_scheduled
level5_unscheduled
level6_unscheduled

shared_hiso_capacity
shared_hiso_used
shared_siso_capacity
shared_siso_used

hiso_mux_unroutable
siso_mux_unroutable
```

row 级调试数据还应记录：

```text
shared_row
source_level
source_local_row
early_stop_hit
hybrid_class
eligibility
class_priority
level_priority_rank
final_action
assigned_core
produced_row
```

这些字段可以直接比较公平、第五级优先和第六级优先三种策略的 BER、资源利用率和两级饥饿情况。

## 16. 验证顺序

建议按以下顺序验证：

1. `Shared HISO=64, Shared SISO=64`，确认共享模式没有因参数、拆分或写回改变基线结果。
2. 固定满配资源，对比 `Fair`、`Level5First`、`Level6First`，理论上三者应得到相同结果。
3. 缩小 SISO 容量，保持 HISO 满配，观察 soft resource sharing。
4. 缩小 HISO 容量，保持 SISO 满配，观察 hard resource sharing。
5. 同时限制 HISO/SISO，比较三种跨级策略的 post-FEC BER。
6. 检查每个 invocation 的资源守恒和 64 行状态守恒。
7. 检查 `EARLY_STOP_BIND_GROUP_SIZE=4` 时不存在跨级 group。
8. 检查 MUX 只改变 core mapping，不改变优先级处理输出的 code 集合；若出现拓扑路由失败，应通过显式补位流程返回优先级模块处理。

## 17. 最终流程定义

第五/六级共享解码的最终定义为：

```text
Level 5 的 32 code + Level 6 的 32 code
  -> 按 [L5 0..31][L6 0..31] 合并，不交错
  -> 64 行分别做 early-stop
  -> bind-group 只在各级内部生效
  -> 未早停行做 hybrid classify-only
  -> 建立资源资格
  -> 跨级统一优先级处理
       - 分类优先级
       - Fair / Level5First / Level6First
       - HISO/SISO 路径选择
       - 预算裁剪和 Unscheduled 决策
  -> HISO MUX 和 SISO MUX 只做 code-to-core 路由
  -> 执行被路由的 HISO/SISO code
  -> 按来源拆回第五级和第六级
  -> 分别 normalize、乘各自 alpha、使用各自 beta
  -> 先写回第五级，再写回第六级
```

这套边界确保“优先级与仲裁”集中在一个模块中，而 MUX 保持纯路由职责，便于后续分别研究跨级优先级算法和具体硬件互连拓扑。
