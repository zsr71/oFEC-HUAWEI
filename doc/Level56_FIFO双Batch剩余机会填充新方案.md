# Level5/6 FIFO 双 Batch 剩余机会填充方案

日期：2026-09-02  
状态：核心调度逻辑已实现并通过回归测试；已完成第一组同种子长帧 BER/FIFO 配对实验，尚需增加独立 seed 验证稳定性。

## 1. 背景

当前 Level5/6 Buffered FIFO 在每个调度时刻最多对一个普通 batch 执行一轮 Group4 调度。一次普通服务内部可以使用多个 Group4 entry，也可以安排多个 HISO/SISO code；“最多一次普通服务”限制的是普通 batch 数量，而不是 code 数量。

当队首 batch 实际用不满 `group4_max_entries=8` 时，当前时刻剩余的 entry 不会交给第二个普通 batch。即使 FIFO 后面仍有可调度工作，这些机会也会闲置。新方案的目标是在不增加单时刻 8-entry 硬件预算的前提下，让 FIFO 中第二个普通 batch 填充第一个普通 batch 用不完的机会。

## 2. 方案正式定义

> 每个调度时刻按照 FIFO 顺序处理普通 batch，普通服务对象最多为两个。第一个普通 batch 优先获得最多 8 次 Group4 entry 机会；若它实际使用不足 8 次，剩余机会可以继续分配给第二个普通 batch。两个 batch 的实际使用次数之和不得超过 8。第二个普通 batch 处理结束后，即使仍有剩余机会，也不再服务第三个普通 batch。

方案本身只有以下三条规则。

### 规则一：每时刻最多服务两个普通 batch

```text
ordinary_batch_count <= 2
```

- 一个调度时刻可以只服务一个普通 batch，也可以服务两个普通 batch。
- 若第一个普通 batch 已用满 8 次 entry，则第二个普通 batch 本时刻不获得服务。
- 即使第二个普通 batch 结束后仍有 entry 剩余，也不继续服务第三个普通 batch。
- FullEarlyStop 不属于普通服务，不占用这两个普通 batch 名额。

### 规则二：严格按照 FIFO 顺序分配机会

第一个普通 batch 的优先级高于第二个普通 batch：

```text
第一个普通 batch 先规划并实际使用 entry
第二个普通 batch 只能使用前者未使用的 entry
```

例如总预算为 8，第一个 batch 可使用 6 次，第二个 batch 可使用 5 次，则实际分配为：

```text
第一个 batch：6
第二个 batch：2
```

不能为了平衡两个 batch 而改成 `4+4`，也不能让第二个 batch 抢占第一个 batch 能够使用的机会。

这里不会越过尚未完成的第一个普通 batch。Group4 调度具有以下性质：第一个普通 batch 如果没有完成，就会用满本时刻全部 8 次机会；只有当它实际使用少于 8 次并已经完成、写回和退休后，新的队首才可能成为第二个普通 batch。

### 规则三：每时刻总 Group4 entry 仍然最多为 8

设第一个、第二个普通 batch 的实际 entry 使用数分别为 `E1` 和 `E2`：

```text
E1 <= 8
E2 <= 8 - E1
E1 + E2 <= 8
```

例如：

| 第一个 batch 可用需求 | 第二个 batch 可用需求 | 实际分配 | 未使用 |
|---:|---:|---|---:|
| 3 | 2 | `3 + 2` | 3 |
| 3 | 7 | `3 + 5` | 0 |
| 8 | 任意 | `8 + 0` | 0 |
| 1 | 7 | `1 + 7` | 0 |

“两个 batch 使用同一个总预算”和“第二个 batch 使用剩余机会”是前三条规则执行后的结果，不作为额外的第四条规则。

## 3. 调度过程

推荐的逻辑过程如下：

```text
remaining_entries = 8
ordinary_batch_count = 0

按 FIFO 顺序检查 batch：
    如果当前 batch 可 FullEarlyStop：
        直接退休
        不占普通 batch 名额
        不消耗 entry
        继续按 FIFO 顺序检查

    如果 ordinary_batch_count == 2：
        结束本时刻普通调度

    用 remaining_entries 作为本 batch 的 entry 上限
    对当前 batch 执行一轮普通规划和解码
    remaining_entries -= 本 batch 实际使用的 entry
    ordinary_batch_count += 1

    如果 batch 完成：
        立即写回全局工作内存
        以 Normal 退休
    否则：
        它必然已经用满本批次预算
        保留在 FIFO 队首
        结束本时刻普通调度

    如果 remaining_entries == 0：
        结束本时刻普通调度

第二个普通 batch 处理后，无论 remaining_entries 是否大于 0：
    本时刻结束，不服务第三个普通 batch
```

## 4. 需要在实现前固定的边界语义

以下内容不改变三条核心规则，但实现和统计时必须保持一致。

### 4.1 什么算一次普通服务

建议只有当 batch 实际使用至少一个 Group4 entry 时，才记作一次普通服务。仅被检查或分类、但没有安排任何 entry，不应增加服务次数。

### 4.2 第二个 batch 的一轮是否允许被剩余预算截断

按照“剩余 5 次机会可以给第二个 batch”的原始定义，第二个 batch 应使用：

```text
min(第二个 batch 当前可用需求, remaining_entries)
```

也就是说，第二个 batch 若可使用 7 次而只剩 5 次，可以在本时刻先使用 5 次，剩余工作留到后续时刻。这是本方案的直接语义，后续实现不应把第二个 batch 限制为“只有完整一轮能装入剩余预算才允许服务”。

### 4.3 FullEarlyStop

沿用现有语义：FullEarlyStop 不占 entry，也不占两个普通 batch 名额。第一普通 batch 完成、写回并退休后，可以继续处理新队首上的一个或多个 FullEarlyStop batch；随后遇到的普通 batch 才作为第二普通 batch。

### 4.4 写回后重新读取

第一普通 batch 完成后，必须先把 Level 5、Level 6 结果写回全局工作内存，再从更新后的全局工作内存读取和分类第二普通 batch。不能在第一普通 batch 写回前缓存第二普通 batch 的输入或分类结果。

实际坐标核验表明，相邻 batch 之间存在 512 个真实重复坐标：前一 batch 的 Level 5 历史写回会成为后一 batch 的 Level 6 当前行输入。因此上述顺序属于实际数据依赖要求。

### 4.5 ForcedEvicted 与窗口边界

现有实现是在一个调度时刻结束时，根据原始队首和 `S_t/window_start` 边界决定是否 ForcedEvicted。新方案中第二个 batch 提前获得服务后，不能简单套用“第二个 batch 也独立触发一次边界淘汰”的逻辑；否则会改变每时刻的缓冲窗口演进规则。

建议实现时保持：

- 一个调度时刻只执行一次全局 `S_t` 更新；
- `completed_batches` 可以是该时刻所有退休 batch 的合计；
- ForcedEvicted 仍以该时刻原始队首和既有边界规则为基准；
- 通过单元测试覆盖“第一 batch 未完成时不服务第二 batch”“两个 batch 均完成”“中间存在 FullEarlyStop”等情况。

### 4.6 HISO/SISO 容量

三条规则直接规定的是 8 个 Group4 entry 总预算。当前 Group4 调度中，entry slot 和 HISO/SISO core 分配存在对应关系。实现第二个 batch 时，应从第一个 batch 已使用的 `entry_slot_offset` 继续编号，并验证同一时刻不会重复占用 core/slot。该约束是三条调度规则落到现有硬件资源模型后的实现结果，不另立为调度规则。

## 5. 旧数据能回答什么

旧结果由单 batch 普通服务策略产生。旧 `level56_schedule_rounds.csv` 能准确给出每次普通服务的 `total_group_entries`，旧 `level56_buffered_times.csv` 能给出当时 FIFO 深度。因此可以直接回答：

1. 原策略有多少普通服务时刻没有用满 8 次 entry；
2. 原策略的平均 entry 利用率和空闲 entry 总量；
3. 没用满时，FIFO 中是否通常还有第二个 batch，亦即新方案是否存在填充对象；
4. 在“每个 batch 的旧工作量和服务结果不因提前执行而改变”的强假设下，静态重排可能达到怎样的 FIFO 收益。

旧数据不能直接回答：

1. 第二个 batch 提前分类和解码后，其真实 pending、早停判定、writeback/history 是否与旧运行相同；
2. 新方案的真实 FIFO 峰值、ForcedEvicted、服务次数分布和 drain 时长；
3. 新方案的 Post-FEC errors/BER。

原因是第二个 batch 被提前服务会改变队列位置、`S_t`、共享 SRAM 读写时序和后续输入。新方案的最终结论必须由新代码在相同 seeds 下配对实测。

## 6. 旧数据中的直接证据：存在大量可填充的 entry 碎片

以下表只使用 `drain=ON` 的完整运行。`输入阶段利用率` 定义为：

```text
输入阶段实际使用 entry 总数 / (输入阶段普通服务时刻数 * 8)
```

“候选填充时刻”要求旧策略在该时刻：

```text
第一个普通 batch 使用 entry < 8
并且 fifo_depth_before >= 2
```

它证明队列中有第二个 batch 可以考虑，但不保证第二个 batch 在新时序下必然能用完所有剩余 entry。

| 旧配置 | 输入阶段普通服务时刻 | 平均 entry/普通服务 | entry 利用率 | 使用不足 8 的时刻 | 候选填充时刻 | 输入阶段空闲 entry |
|---|---:|---:|---:|---:|---:|---:|
| 3.065 dB，长帧，R_buf=16000，0/8 | 16851 | 5.482 | 68.52% | 65.65% | 65.30% | 42436 |
| 3.065 dB，长帧，R_buf=16000，8/8 | 16851 | 4.959 | 61.98% | 81.89% | 81.35% | 51252 |
| 3.09 dB，短帧，R_buf=1600，0/8 | 8404 | 4.570 | 57.12% | 85.52% | 84.72% | 28828 |
| 3.06 dB，短帧，R_buf=1600，0/8 | 8404 | 6.322 | 79.02% | 48.63% | 47.86% | 14104 |
| 3.06 dB，短帧，R_buf=3200，0/8 | 8404 | 6.137 | 76.71% | 51.45% | 50.68% | 15657 |
| 3.065 dB，长帧，R_buf=6400，0/8 | 16851 | 5.704 | 71.30% | 61.76% | 61.41% | 38696 |

直接结论：

- 所有已保留的 `drain=ON` 配置都有明显的 entry 碎片；输入阶段 entry 利用率约为 57%～79%。
- 在约 48%～85% 的输入阶段普通服务时刻，第一个 batch 没用满 8 次且 FIFO 中至少还有第二个 batch。
- 因而新方案针对的是旧运行中实际存在、而且占比很高的空闲机会，不是只在少数尾部时刻发生的特殊情况。
- 3.06 dB 下单次工作更重，旧策略已经使用约 77%～79% 的 entry 能力，预计填充收益空间小于 3.09 dB 和 3.065 dB 的轻负载样本，但仍不是零。

### 6.1 五组独立 seed 的稳定性

在固定 `Eb/N0=3.065 dB、60,014,592 input bits、R_buf=16000、drain=ON、Group4 entries=8` 的五组独立 seed 中：

| 模式 | 输入阶段 entry 利用率，五 seed 范围 | 使用不足 8 的时刻，五 seed 范围 | 候选填充时刻，五 seed 范围 |
|---|---:|---:|---:|
| 0/8 | 67.44%～68.52%，平均 68.28% | 65.01%～68.05%，平均 66.02% | 64.70%～68.01%，平均 65.79% |
| 8/8 | 60.90%～62.60%，平均 61.95% | 81.03%～84.40%，平均 82.50% | 79.89%～83.89%，平均 81.84% |

五组 seed 的范围很窄，说明 entry 碎片不是单个随机帧的偶然现象。特别是 8/8 虽然旧 FIFO 峰值低于 0/8，但按 entry 统计的碎片更多，因此双 batch 填充仍有明显的资源利用率提升空间。

## 7. 离线重排预测：可能的 FIFO 收益

为得到一个数量级判断，可以对旧数据做静态工作量重排：

1. 保留每个 batch 在旧运行中的普通服务 entry 序列和总工作量；
2. 每个时刻到达关系不变；
3. FIFO 前一 batch 优先；
4. 第一 batch 使用后，第二 batch 最多使用剩余 entry；
5. 每时刻最多两个普通 batch、合计最多 8 次 entry；
6. 第二 batch 的旧工作块允许按剩余预算截断；
7. 假设提前服务不改变早停、解码结果、writeback/history 和未来工作量；
8. 不重新模拟真实 `S_t`/ForcedEvicted/共享 SRAM 状态。

第 7、8 条是很强的假设，因此以下数据是“静态工作量重排预测”，不是新方案仿真结果。

### 7.1 五 seed 的预测范围

| 模式 | 旧 FIFO 峰值 | 静态预测峰值 | 旧输入结束积压 | 静态预测输入结束积压 | 旧 drain 时刻 | 静态预测 drain 时刻 |
|---|---:|---:|---:|---:|---:|---:|
| 0/8，五 seed 范围 | 4037～4590 | 333～730 | 4037～4589 | 176～693 | 5460～6420 | 191～684 |
| 0/8，五 seed 平均 | 4415.8 | 529.0 | 4415.2 | 324.0 | 5907.2 | 344.4 |
| 8/8，五 seed 范围 | 1552～2027 | 104～141 | 1551～2026 | 0～35 | 1799～2342 | 0～28 |
| 8/8，五 seed 平均 | 1804.2 | 121.8 | 1803.4 | 12.8 | 2041.4 | 15.2 |

按五 seed 平均，静态模型预测：

- 0/8：FIFO 峰值约下降 88.0%，输入结束积压约下降 92.7%，drain 时刻约下降 94.2%；
- 8/8：FIFO 峰值约下降 93.2%，输入结束积压和 drain 在静态模型中接近消失。

这些幅度应视为乐观预测。它们的价值是说明：旧运行中的总 entry 工作量从容量角度看并非必然造成数千 batch 的积压，原策略的“单时刻只允许一个普通 batch”确实产生了严重的容量碎片。它们不能代替新实现实测。

### 7.2 为什么静态预测可能过于乐观

- 第二个 batch 提前读取共享 SRAM 后，输入可能不同于旧运行到达队首时的输入；
- 提前 writeback 会改变历史和后续 batch 的状态；
- batch 的下一次可调度需求可能不是旧 entry 序列的简单平移；
- FullEarlyStop、Normal、ForcedEvicted 的退休原因可能改变；
- 一个时刻完成两个 batch 会改变 `S_t` 更新，新窗口轨迹不再等同旧运行；
- 真实 Group4 规划的第二 batch 可能受剩余 slot/core 组合限制，未必能把所有数字意义上的余量填满。

因此，不应在新代码运行前使用上述预测值确定最终 FIFO 深度，也不应据此预测 BER。

## 8. 对新方案收益的当前判断

### 高可信结论

1. 新方案具有真实且稳定的填充机会。多 seed 旧数据中，0/8 约 66%、8/8 约 82.5% 的输入阶段普通服务时刻没有用满 8 次。
2. 大多数这些时刻 FIFO 中确实还有第二个 batch：0/8 候选填充时刻平均约 65.8%，8/8 平均约 81.8%。
3. 新方案不增加每时刻 8-entry 峰值能力，主要提升平均利用率并缓解由 batch 边界造成的资源碎片。
4. 预期方向是降低 FIFO 增长率、峰值、输入结束积压和 drain 时长；在小缓冲配置下还可能减少 ForcedEvicted。

### 中等可信判断

1. 3.065 dB 长帧的收益可能很大，因为旧输入阶段长期积压，且候选填充时刻占比高。
2. 8/8 的 FIFO 原本较浅，但 entry 碎片比例比 0/8 更高，因此仍可能明显受益。
3. 3.06 dB 的单 batch 工作量更重，旧 entry 利用率更高，因此相对收益预计小于 3.09/3.065 dB 的轻负载样本。

### 目前不能声称的结论

1. 不能声称新方案一定把 FIFO 峰值降到静态预测值；
2. 不能声称新方案一定消除连续输入下的正漂移；
3. 不能声称 Post-FEC BER 不变、改善或恶化；
4. 不能用旧 `drain=OFF` 的 BER 推断新方案 BER。

## 9. 实现后的验证方案

建议保留旧策略作为可切换基线，例如：

```text
LEVEL56_FIFO_MAX_ORDINARY_BATCHES_PER_TIME = 1  // 旧策略
LEVEL56_FIFO_MAX_ORDINARY_BATCHES_PER_TIME = 2  // 新策略
```

首轮配对实验固定：

```text
FIFO = ON
Group4 = ON
group4_max_entries = 8
drain = ON
相同 Eb/N0、R_buf、帧长、bit seed、channel seed
仅改变每时刻最大普通 batch 数：1 vs 2
```

至少覆盖：

| 优先级 | 配置 | 目的 |
|---:|---|---|
| 1 | 3.065 dB，长帧，R_buf=16000，0/8，五 seeds | 验证 FIFO 峰值、漂移、drain 与 BER |
| 2 | 3.065 dB，长帧，R_buf=16000，8/8，五 seeds | 验证 8/8 的高碎片是否转化为真实收益 |
| 3 | 3.09 dB，R_buf=1600，0/8，drain=ON | 验证原本已无 ForcedEvicted 的轻负载性能 |
| 4 | 3.06 dB，R_buf=1600/3200，0/8，drain=ON | 验证重负载下是否减少 ForcedEvicted |

新日志至少需要记录：

```text
ordinary_batch_count
ordinary_batch_1 / ordinary_batch_2
entries_used_batch_1 / entries_used_batch_2
entries_used_total / entries_unused
第二 batch 是否被剩余预算截断
每个 batch 的服务次数和退休原因
fifo_depth_before / fifo_depth_after
S_t / S_next
ForcedEvicted / censored
Pre-FEC 和 Post-FEC errors、BER、比较总比特数
```

验收条件：

- 每时刻普通 batch 数不超过 2；
- 每时刻两个普通 batch 的 entry 合计不超过 8；
- 第一 batch 未用完之前，第二 batch 不占用 entry；
- 第二 batch 后不服务第三个普通 batch；
- entry slot/core 不重复；
- `drain=ON` 的完整 BER 运行满足 `censored=0`；
- 对无容量淘汰配置，重点检查 `ForcedEvicted=0`；
- 用同 seed 的旧/新方案逐项比较 FIFO 峰值、输入结束积压、drain、服务次数、退休原因和 BER。

## 10. 当前结论

旧数据足以证明新方案具有明确的资源利用率动机：原策略在大量调度时刻只使用了 8 个 entry 中的一部分，而 FIFO 中通常还有第二个 batch。静态工作量重排显示 FIFO 峰值可能从“数千 batch”下降到“数百甚至百级”，但这一数字依赖于状态不变的强假设，只能作为新方案值得实现和测试的依据。

实现后的短帧冒烟测试使用 `3.065 dB、710400 bits、R_buf=32、0/8、Group4 entries=8、drain=ON`。共观察到 34 个时刻服务两个普通 batch，其中输入阶段 21 个、帧尾 drain 阶段 13 个；每个时刻实际使用 entry 总数均不超过 8，第二 batch 的 slot 偏移等于第一 batch 的实际使用数，FIFO 最终排空。该运行用于验证控制流和日志，不用于得出正式 BER/FIFO 收益结论。

### 10.1 第一组同种子长帧配对结果

第一组正式配对固定以下配置，仅改变普通调度策略：

```text
Eb/N0 = 3.065 dB
信息比特数 = 60014592
bit seed = 20260319
channel seed = 3182026
FIFO = ON
R_buf = 16000
drain = ON
HISO/SISO = 0/8
Group4 entries = 8
```

| 指标 | 原单 Batch 策略 | 双 Batch 策略 | 变化 |
|---|---:|---:|---:|
| Pre-FEC errors / bits | 1292037 / 58608000 | 1292037 / 58608000 | 完全一致 |
| Post-FEC errors / bits | 11 / 58608000 | 11 / 58608000 | 完全一致 |
| Post-FEC BER | 1.87688e-7 | 1.87688e-7 | 完全一致 |
| Level56 entry 总工作量 | 125500 | 125500 | 完全一致 |
| FIFO 峰值 | 4476 | 360 | 下降 91.96% |
| 输入结束积压 | 4475 | 197 | 下降 95.60% |
| drain 时刻数 | 6074 | 202 | 下降 96.67% |
| 输入加 drain 总时刻数 | 22971 | 17099 | 下降 25.56% |
| 最小 `S_t` | 3852 | 15282 | 提高 11430 行 |
| ForcedEvicted | 0 | 0 | 均为 0 |
| Normal / FullEarlyStop | 16781 / 116 | 16781 / 116 | 完全一致 |

在只统计“输入阶段发生普通服务”的时刻时，entry 利用率由 68.52% 提高到 92.05%，提高 23.53 个百分点。双 Batch 运行中共有 11125 个时刻服务了两个普通 batch，其中输入阶段 10984 个、drain 阶段 141 个；3012 个双 Batch 时刻中两个普通 batch 都完成。

对全部 11125 个双普通 Batch 时刻执行轨迹约束检查，以下违例数量均为 0：

- 第一普通 batch 未完成却继续服务第二普通 batch；
- 第二普通 batch 的 slot 偏移不等于第一普通 batch 的实际使用数；
- 第二普通 batch 的预算不等于 `8 - 第一普通 batch 实际使用数`；
- 两个普通 batch 的使用数之和与时刻总使用数不一致；
- 一个时刻的总使用数超过 8。

这组结果说明，对该 seed，双 Batch 策略把相同的 125500 次 Group4 entry 工作更紧凑地安排到了输入阶段，显著减少了积压和 drain，同时没有改变 Pre-FEC/Post-FEC 错误数。新旧日志中展示的前 10 个 Post-FEC 错误位置也一致。由于目前只有一个正式长帧 seed，不能据此宣称所有 seed 的 BER 都必然不变；下一步仍应完成其余独立 seed 配对。

最终应以新实现、相同 seed、`drain=ON` 的配对结果回答三个问题：

1. entry 利用率实际提高了多少；
2. FIFO 峰值、输入结束积压、drain 和 ForcedEvicted 实际下降了多少；
3. 提前服务第二个 batch 是否改变最终 BER。

## 附录 A：可复查数据与脚本

本页统计脚本：

```text
analysis_fifo_two_batch_potential.py
```

运行方式：

```bash
python3 analysis_fifo_two_batch_potential.py
```

脚本的第一张输出表是旧轨迹的直接统计；第二张输出表带有 `Optimistic replay; not a decoder simulation` 标题，是乐观静态重排，不能当作新方案仿真结果。

主要单次配置的数据源位于 `data/level56_schedule/`，包括：

```text
ofec_single_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso0_siso8_schedgroup4_level56_schedule_rounds.csv
ofec_single_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso8_siso8_schedgroup4_level56_schedule_rounds.csv
ofec_single_ebn03.09_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_schedule_rounds.csv
ofec_single_ebn03.06_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_schedule_rounds.csv
ofec_single_ebn03.06_fifo1_rbuf3200_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_schedule_rounds.csv
ofec_single_ebn03.065_nbits60014592_fifo1_rbuf6400_drain1_hiso0_siso8_schedgroup4_level56_schedule_rounds.csv
```

每个 `schedule_rounds.csv` 通过同 stem 的 `level56_buffered_times.csv` 对齐输入阶段、FIFO 深度和旧峰值。五 seed 统计使用文件名中含有以下字段的成对数据：

```text
ebn03.065_nbits60014592_bitseed*_chseed*_fifo1_rbuf16000_drain1_hiso0_siso8_g4entries8_schedgroup4
ebn03.065_nbits60014592_bitseed*_chseed*_fifo1_rbuf16000_drain1_hiso8_siso8_g4entries8_schedgroup4
```
