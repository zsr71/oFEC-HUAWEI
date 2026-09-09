# Level56 FIFO 多 Batch 连续服务与译码延迟性能验证方案

日期：2026-09-09
状态：方案定义；尚未创建实验 worktree，尚未实现，尚未产生 BER/FIFO 实测结果。

## 1. 目标

本方案以 `feature/level56-buffered-fifo-scheduling` 为基线，在 Level 5/6
共享 FIFO 中加入细时钟级的译码延迟约束，并验证它对下列结果的影响：

- Post-FEC BER 与错误位置；
- FIFO 深度、pending 数量和帧尾未完成 batch 数量；
- batch 完成、Full-EarlyStop 退休和 ForcedEvicted 数量；
- 一个到达时刻内 8 次 Group4 机会的利用情况；
- 延迟阻塞造成的空闲时钟数量。

本阶段追求的是延迟约束对最终 BER 和 FIFO 行为的有效反馈，不追求逐寄存器、
逐端口复刻真实硬件流水线。SISO/HISO 数值译码仍在被调度时立即计算并写回；
本阶段不引入“在途译码结果队列”，也不模拟结果在若干时钟后才写回。

## 2. 基线与实验边界

方案基线为：

```text
仓库：/home/zsr71/projects/newcode_level56_buffered_fifo
分支：feature/level56-buffered-fifo-scheduling
基线提交（编写本文档时）：c1d782e7c8023538d7c1c30228a918a98459a74b
```

当前基线的一个 Level 5/6 FIFO batch 是完整的 64-code batch：

```text
Level 5：32 code = 8 个四码字组
Level 6：32 code = 8 个四码字组
合计   ：64 code = 16 个四码字组
```

本方案只修改 Level 5/6 buffered FIFO 路径。FIFO 关闭、Level 5/6 shared 关闭、
temporal lookahead 等其它路径必须保持原行为。

## 3. 时间和 Group4 机会的统一口径

### 3.1 一个到达时刻包含 8 个细时钟

每隔 8 个细时钟到达一个新的 64-code batch：

```text
B1 到达：绝对时钟 0
B2 到达：绝对时钟 8
B3 到达：绝对时钟 16
B4 到达：绝对时钟 24
B5 到达：绝对时钟 32
```

到达时刻 `t` 对应的绝对时钟区间为：

```text
[8t, 8t+8)
```

区间内有 8 个时钟槽位：

```text
槽位 0、1、2、3、4、5、6、7
```

### 3.2 每个细时钟只对应一个四码字组机会

每个细时钟只能使用一次 Group4 机会，即只处理一个四码字组。该组或者来自
Level 5，或者来自 Level 6，不会在同一时钟同时读取 Level 5 和 Level 6 两组。

Level 5/6 合计 16 个四码字组，共同竞争一个到达时刻内的 8 次 Group4 机会。
当前 Group4 负载排序、多轮选择、组内 SISO/HISO 动作选择继续沿用现有调度器。

延迟模型只关心某个 batch 实际占用了哪些 Group4 时钟槽位，不区分该槽位中：

- 只安排了 SISO；
- 只安排了 HISO；
- 同时安排了 SISO 和 HISO。

只要一个 Group4 entry 被实际使用，对延迟模型而言就消耗一个细时钟。

### 3.3 Group4 机会槽位与物理组编号不同

槽位 `0..7` 是当前到达时刻内的时间位置，不是 16 个候选四码字组的物理编号。
例如：

```text
槽位 0 可能处理 Level 5 的组 3
槽位 1 可能处理 Level 6 的组 1
槽位 2 可能再次处理仍有 pending code 的某个组
```

延迟起点使用的是最后占用的时间槽位，而不是被选中组的物理编号。

## 4. 一个时刻可以连续服务多个 batch

### 4.1 取消固定普通 batch 数量上限

本方案不再规定一个到达时刻最多服务一个或两个普通 batch。一个时刻内可以按
FIFO 优先级连续服务任意多个已经到达的普通 batch。

实际可服务的 batch 数由以下条件自然限制：

- 当前区间总共只有 8 个 Group4 时钟；
- 每个 batch 实际消耗的 Group4 时钟数；
- 当前 batch 是否受到译码延迟保护；
- FIFO 中是否还有已经到达并可处理的 batch。

如果每个 batch 都只使用一个 Group4 时钟，一个到达时刻最多自然服务 8 个
普通 batch；这不是额外的 `max_batches=8` 规则，而是 8 个总时钟推导出的上限。

### 4.2 剩余时钟继续交给后续 batch

假设某个区间开始时 FIFO 中已有 B1、B2、B3，B1 使用槽位 `0..5`：

```text
槽位 0：B1
槽位 1：B1
槽位 2：B1
槽位 3：B1
槽位 4：B1
槽位 5：B1
槽位 6：若 B2 已到达且不受保护，则服务 B2
槽位 7：继续服务 B2，或在 B2 本轮结束后服务下一个可处理 batch
```

只有当 B2 尚未到达、FIFO 没有其它可处理 batch，或者下一个按顺序应处理的
batch 受到延迟保护时，剩余槽位才会空闲。

### 4.3 同一时刻的临时前移不改变下一时刻的 FIFO 优先级

一个 batch 本时刻获得一轮服务但尚未完成时，允许本时刻的服务游标继续访问
它后面的 batch，从而使用剩余 Group4 时钟。下一到达时刻仍从最早的未完成
batch 重新开始。

例如：

```text
时刻 t：
  B1 使用槽位 0..2，仍未完成
  B2 使用槽位 3..4，仍未完成
  B3 使用槽位 5，完成并退休
  B4 使用槽位 6..7

时刻 t+1：
  仍先从最早未完成的 B1 开始，再到 B2、B4……
```

因此需要区分：

- FIFO 中永久保持的到达与优先级顺序；
- 当前 8 时钟区间内临时向后推进的服务游标。

## 5. 译码延迟起点与解除时刻

### 5.1 延迟参数

第一阶段建议增加：

```cpp
std::size_t LEVEL56_SISO_DECODER_LATENCY = 0;
```

单位为细时钟。性能对齐实验的主要值为：

```text
LEVEL56_SISO_DECODER_LATENCY = 16
```

虽然沿用 `SISO_DECODER_LATENCY` 名称，但本方案的延迟起点由最后使用的 Group4
机会决定，不检查该机会内部最终执行的是 SISO 还是 HISO。文档、日志和测试必须
明确这一实验语义，避免把它误解为严格的 SISO 结果返回模型。

### 5.2 起点取 batch 最近一次实际使用的最后槽位

如果当前到达时刻的绝对起始时钟为 `T`，某 batch 本轮最后使用的槽位为 `s`，
则该 batch 最近一次 Group4 处理时钟为：

```text
last_used_cycle = T + s
```

延迟解除时钟定义为：

```text
release_cycle = last_used_cycle + LEVEL56_SISO_DECODER_LATENCY
```

不额外加 1。

判断规则为：

```text
current_cycle < release_cycle：仍受保护
current_cycle >= release_cycle：保护已经解除，当前时钟即可使用
```

这是本方案固定采用的时钟编号口径。例如 B1 使用槽位 `0..5`：

```text
B1 最后使用时钟 = 5
延迟起点         = 5
延迟             = 16
解除时钟         = 5 + 16 = 21
时钟 21 起允许受影响的后续 batch 处理
```

此例中 B4 在时钟 24 才到达，所以 B1 的这次保护不会让 B4 实际等待。

### 5.3 多轮服务时更新起点

一个 batch 可能跨多个到达时刻接受服务。每次实际使用 Group4 机会后，立即用本轮
最后使用的绝对时钟更新该 batch 的：

```text
last_used_cycle
release_cycle
```

例如：

```text
B1 第一轮最后使用时钟：7
B1 第二轮最后使用时钟：11

最终采用：
last_used_cycle = 11
release_cycle   = 11 + 16 = 27
```

如果本轮没有实际使用 Group4 entry，则不能更新延迟起点。Full-EarlyStop 零时钟
退休、仅分类、仅检查、受保护而未服务等情况都不建立新的延迟。

### 5.4 延迟状态立即生效，不等待 batch 退休

同一 8 时钟区间可以继续服务后续多个 batch，因此源 batch 每次使用完本轮
Group4 机会后，必须立即更新自己的延迟状态。不能等源 batch 完成或退休以后才建立
保护。

例如：

```text
B1 使用槽位 0..2，但仍未完成
B2 使用槽位 3
B3 使用槽位 4
接下来准备处理 B4
```

B4 已经是 B1 后面的第 3 个 batch，所以在检查 B4 时必须立即看到 B1 在槽位 2
建立的延迟约束。

## 6. 哪些后续 batch 受保护

沿用参考实验的抽象依赖距离：源 batch 后面的第 1、2 个 batch 可以继续；第 3 个
及更远 batch 在保护解除前不能处理。

对于源 batch `Bi` 和候选 batch `Bj`：

```text
j <= i + 2：不受 Bi 这条保护限制
j >= i + 3：受 Bi 这条保护限制
```

例如：

```text
B1 建立保护：
  B2 可继续
  B3 可继续
  B4 及以后需要检查 B1 的 release_cycle
```

延迟约束只限制普通 Group4 服务。全 EarlyStop batch 允许沿用现有零普通预算的完成
与退休逻辑，但不得借此让普通未完成 batch 越过仍有效的保护。

## 7. 多条保护同时存在

候选 batch 可能同时受到多个更早 batch 的保护。应找出所有满足依赖距离条件的
源 batch，并取最晚解除时钟：

```text
effective_release_cycle = max(所有适用保护的 release_cycle)
```

例如：

```text
B1 解除时钟 = 28
B2 解除时钟 = 31
B3 解除时钟 = 35
当前候选 = B5

B5 相对 B1 是第 4 个后续 batch：受 B1 约束
B5 相对 B2 是第 3 个后续 batch：受 B2 约束
B5 相对 B3 是第 2 个后续 batch：不受 B3 约束

B5 最早可处理时钟 = max(28,31) = 31
```

## 8. 延迟在区间中间解除

延迟保护不能以整个 `service_time` 为粒度粗略跳过。若保护在当前 8 时钟区间中间
解除，应立即使用区间内剩余的 Group4 机会。

例如 B4 在时钟 24 到达，有效解除时钟为 28：

```text
时钟 24：受保护
时钟 25：受保护
时钟 26：受保护
时钟 27：受保护
时钟 28：保护解除，可使用槽位 4
时钟 29：可使用槽位 5
时钟 30：可使用槽位 6
时钟 31：可使用槽位 7
```

因此 B4 在该区间仍有 4 次 Group4 机会，不能把整个 `[24,32)` 直接判为空闲。

实现上可维护绝对时钟游标：

```text
interval_begin = service_time * 8
interval_end   = interval_begin + 8
cursor         = interval_begin
```

如果候选 batch 的 `effective_release_cycle > cursor`，先把 `cursor` 推进到解除时钟；
只要仍满足 `cursor < interval_end`，就用剩余槽位继续调度。

## 9. 一个到达时刻内的建议调度流程

```text
1. 在 interval_begin 到达一个新的 64-code batch，并进入 FIFO。

2. 令 cursor = interval_begin，从最早未完成 batch 开始设置本时刻服务游标。

3. Full-EarlyStop batch 按现有规则零普通时钟完成、写回和退休，然后继续检查后续
   batch。

4. 对当前普通 batch 计算 effective_release_cycle。

5. 若 effective_release_cycle > cursor，则把 cursor 推进到解除时钟；被跳过的差值
   计为 latency_blocked_cycles。

6. 若 cursor >= interval_end，则本到达时刻结束。

7. 剩余 Group4 预算为 interval_end - cursor；当前调度的 entry_slot_offset 为
   cursor - interval_begin。

8. 对当前 batch 运行一轮现有 Group4 调度，最大 entry 数不得超过剩余预算。

9. 根据实际使用的 Group4 entry 数 E 消耗 E 个连续时钟，并立即更新该 batch 的
   last_used_cycle 和 release_cycle。

10. 如果 batch 完成则退休；如果未完成则保留在原 FIFO 位置，但本时刻服务游标可
    临时继续访问后面的 batch。

11. 如果 E=0 且不是 Full-EarlyStop，则停止或显式记录 no-progress，防止死循环。

12. 重复步骤 3～11，直到 8 个时钟用完、没有已到达 batch，或者所有剩余候选在
    本区间结束前均无法处理。

13. 下一个到达时刻重新从最早未完成 batch 开始。
```

当前 Group4 调度器已经具有 `max_group_entries` 和 `entry_slot_offset` 参数，可以用于
表达“延迟在中间解除以后，只使用后半段槽位”的行为。新增逻辑不应重新实现组内
分类、优先级、HISO/SISO 选路或写回。

## 10. 关键例子

### 10.1 B1 使用槽位 0..5，后续批次已在 FIFO 中

假设当前区间开始时 B1、B2、B3 已经到达，且均不受旧保护：

```text
槽位 0..5：B1 使用 6 次 Group4 机会
槽位 6..7：继续服务 B2
```

B1 的延迟状态为：

```text
last_used_cycle = 5
release_cycle   = 5 + 16 = 21
```

B2 必须按照它自己最后实际使用的槽位更新独立的延迟状态。

### 10.2 同一时刻连续服务 B1、B2、B3，B4 被延迟挡住

假设区间为 `[100,108)`，FIFO 中已有 B1～B6：

```text
时钟 100..102：B1 使用 3 次机会
时钟 103     ：B2 使用 1 次机会
时钟 104     ：B3 使用 1 次机会
```

B1 最后使用时钟为 102：

```text
B1.release_cycle = 102 + 16 = 118
```

B4 是 B1 后面的第 3 个 batch。准备在时钟 105 服务 B4 时：

```text
105 < 118
```

所以 B4 不能处理。B5、B6 不得越过 B4。当前区间结束于 108，早于 118，因此
时钟 105..107 计为延迟阻塞，本时刻结束。

停止原因不是“已经服务了 3 个 batch”，而是 B4 受到 B1 的延迟保护。

### 10.3 保护在区间中间解除后继续服务多批次

假设区间为 `[24,32)`，B4 的有效解除时钟为 28：

```text
时钟 24..27：延迟阻塞
时钟 28..29：B4 使用 2 次机会
时钟 30     ：B5 使用 1 次机会
时钟 31     ：B6 使用 1 次机会
```

本时刻可以在等待结束后连续服务 B4、B5、B6，不设置普通 batch 数量上限。

### 10.4 延迟为 0 的边界

当：

```text
LEVEL56_SISO_DECODER_LATENCY = 0
```

则：

```text
release_cycle = last_used_cycle
```

由于服务游标在使用 entry 后已经向后推进，任何后续 batch 都不应因为零延迟新增
空闲时钟。新分支在延迟为 0 时必须与“多 batch 连续填满 8 个 Group4 时钟、但不加
延迟保护”的对应基线完全一致。

## 11. FIFO 边界和帧尾 drain

### 11.1 FIFO 边界

延迟会减少有效 Group4 时钟，可能使 pending 和 FIFO 深度增长。触达既有缓冲边界
时，ForcedEvicted 仍按当前 FIFO 方案的物理窗口边界和退休语义执行；延迟实现不能
悄悄改变 forced eviction 的数据定义。

### 11.2 帧尾 drain

需要分别验证：

```text
LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END = 0
LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END = 1
```

drain 阶段不再引入新 batch，但细时钟和延迟保护仍应继续前进。保护期内可以产生
空闲时钟；保护解除后继续服务 FIFO，直到所有可完成 batch 退休。若实现中保留独立
延迟记录，drain 结束时也应确认不存在会影响未完成 batch 的遗留保护状态。

## 12. 参数与运行入口建议

建议新增参数：

```cpp
std::size_t LEVEL56_SISO_DECODER_LATENCY = 0;
```

建议在 `ofec_single` 中增加环境变量覆盖：

```text
LEVEL56_SISO_DECODER_LATENCY=0
LEVEL56_SISO_DECODER_LATENCY=8
LEVEL56_SISO_DECODER_LATENCY=16
LEVEL56_SISO_DECODER_LATENCY=24
LEVEL56_SISO_DECODER_LATENCY=32
```

运行标签和输出文件名必须包含延迟值，例如：

```text
_sisolat0
_sisolat16
```

避免不同延迟实验覆盖同名 CSV 或 BER 输出。

## 13. 观测字段

现有 `level56_buffered_times.csv` 建议增加：

```text
interval_begin_cycle
interval_end_cycle
first_used_slot
last_used_slot
group_entries_used
ordinary_batches_served
latency_blocked_cycles
cumulative_latency_blocked_cycles
candidate_batch
effective_release_cycle
blocking_source_batches
fifo_depth_after
pending_after
```

逐 batch 建议记录：

```text
batch_id
arrival_cycle
service_attempts
group_entry_cycles_used
last_used_cycle
release_cycle
completion_cycle
retirement_reason
```

这些字段只用于观测，不得反向参与调度决策。

## 14. 实现分层建议

建议按以下提交拆分：

1. 增加参数、细时钟字段、batch 延迟字段和 CSV 字段；默认延迟为 0。
2. 将单普通 batch 服务改成同一时刻多 batch 连续服务，不设置固定批次数上限。
3. 加入第 3 个后续 batch 开始受限的细时钟延迟保护，以及区间中途解除逻辑。
4. 增加确定性回归场景并完成延迟 0 对照。
5. 运行固定配置和多 seed BER/FIFO 性能实验，记录结论。

本阶段不实现：

- SISO 结果在途队列；
- 延迟结束后才写回 LLR/history；
- SISO 核持续占用模型；
- SRAM 多写口返回冲突；
- 逐寄存器传输级硬件流水线。

## 15. 必须覆盖的确定性测试

至少覆盖：

1. 一个 batch 恰好使用 8 个槽位。
2. B1 使用 `0..5` 后，已到达且不受保护的 B2 使用 `6..7`。
3. B1、B2、B3 每批只使用少量槽位，同一时刻连续服务三批以上。
4. 不允许因为固定的普通 batch 数量上限提前停止。
5. B1 最后使用槽位 5、延迟 16，解除时钟严格为 21。
6. 第 1、2 个后续 batch 不受源 batch 保护。
7. 第 3 个后续 batch 在解除前受保护，在解除时钟立即恢复。
8. 保护在 8 时钟区间中间解除，只消耗前半段，后半段继续服务。
9. 多条适用保护取最晚解除时钟。
10. batch 多轮服务后使用最近一次最后槽位更新延迟起点。
11. Full-EarlyStop 不消耗普通 Group4 时钟，也不建立新延迟。
12. 延迟阻塞时后续普通 batch 不得越过当前受保护 batch。
13. 延迟为 0 时不引入额外空闲时钟。
14. no-progress 场景不会死循环。
15. drain 阶段没有新 batch 到达，但会推进细时钟并最终清空 FIFO。

## 16. BER 性能验证矩阵

第一轮固定相同 Eb/N0、帧长、比特种子、信道种子、量化参数、早停参数、
HISO/SISO 配置、Group4 调度、缓冲深度，只改变延迟：

| 组别 | 多 batch 连续服务 | 延迟/细时钟 | drain | 用途 |
|---|---|---:|---|---|
| A | 开 | 0 | 关 | 在线无延迟基线 |
| B | 开 | 16 | 关 | 在线延迟效果 |
| C | 开 | 0 | 开 | 所有 batch 完成的无延迟参考 |
| D | 开 | 16 | 开 | 排除帧尾未完成后的延迟对照 |
| E | 开 | 8/24/32 | 关、开 | 延迟敏感性扫描 |

至少比较：

- Pre-FEC、Quantized-hard 和 Post-FEC BER；
- Post-FEC 错误位置文件；
- completed、Full-EarlyStop、ForcedEvicted；
- 最大/平均 FIFO 深度；
- 输入阶段和帧尾 pending batch 数量；
- Group4 entry 利用率；
- 延迟阻塞细时钟数量和比例；
- batch 到达至完成的平均值、最大值和分位数。

第一轮固定样本通过后，再在多个 seed 和多个 Eb/N0 点重复。延迟效果与早停率及
Group4 负载直接相关，不能用单个 seed、单个信噪比作为普遍结论。

## 17. 验收标准

### 17.1 结构验收

- 每个到达时刻严格只有 8 个 Group4 细时钟；
- 每个细时钟只对应一个四码字组机会；
- 同一时刻普通 batch 数量不设固定上限；
- 已使用槽位和延迟阻塞槽位之和不超过 8；
- 下一时刻重新从最早未完成 batch 开始；
- 不发生后续普通 batch 越过受保护 batch；
- 延迟 16 的起点和解除时钟符合第 5 节定义。

### 17.2 回归验收

- FIFO/延迟相关功能关闭时，与原分支行为完全一致；
- 新多 batch 连续服务打开但延迟为 0 时，有独立且可复现的无延迟基线；
- 观测开关关闭和打开时，BER、错误位置、调度动作与 FIFO 生命周期一致；
- 既有 Level 5/6 shared regression check 全部通过。

### 17.3 性能结论验收

最终报告必须把以下影响分开：

- 多 batch 连续使用剩余 Group4 机会带来的吞吐改善；
- 延迟保护消耗 Group4 时钟带来的吞吐损失；
- FIFO 边界和 ForcedEvicted 对 BER 的影响；
- 帧尾未 drain 与完整 drain 对最终 BER 的影响。

只有同时给出 BER、FIFO、pending、Group4 利用率和延迟阻塞时钟，才对实验结果
作出因果解释。

## 18. 当前正式方案摘要

```text
一个 FIFO batch = Level 5 的 32 code + Level 6 的 32 code = 64 code。
64 code 按每组 4 code 分为 16 个组。
每隔 8 个细时钟到达一个新 batch。
每个细时钟只能处理一个四码字组。
一个到达时刻可以连续服务任意多个已经到达的普通 batch。
不设置最多一个、两个或其它固定普通 batch 数量限制。
一个 batch 用完部分槽位后，剩余槽位继续给后续可处理 batch。
下一到达时刻重新从最早未完成 batch 开始。
延迟起点取 batch 最近一次实际使用的最后 Group4 时钟。
不区分该组内执行的是 SISO 还是 HISO。
解除时钟 = 最后使用时钟 + SISO_DECODER_LATENCY；不额外加 1。
当前目标延迟为 16 个细时钟。
源 batch 后面的第 1、2 个 batch 可继续；第 3 个及以后受保护。
保护在区间中途解除时，立即使用剩余时钟继续服务。
译码和写回仍立即完成；本阶段不实现延迟写回。
最终通过 BER、FIFO 深度、pending、forced eviction 和阻塞时钟共同验证效果。
```
