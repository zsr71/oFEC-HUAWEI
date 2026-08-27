# Level 5/6 Buffered FIFO 与独立充足资源基线的等价性验证计划

## 1. 目的与待证命题

本计划验证的不是“缓冲区是否降低了 BER”这一宽泛问题，而是两个彼此独立的命题。

### 1.1 命题 5.1：译码机会覆盖性

在 Buffered FIFO 方案中，若缓冲区足够长，且一个 batch 未因物理窗口边界而被强制顶出，则：

- 每个 EarlyStop code 都执行原有 `EarlyStopAction` 并完成写回；
- 每个初始为非 EarlyStop 的 code 最终至少获得一次 `HisoDecode` 或 `SisoDecode`；
- 已完成 code 不会重复占用资源或重复写回；
- 不存在最终 `pending`、`forced_evicted`、`ActionFailed` 或未完成写回的 code。

该命题回答：新方案是否真正消除了“某些待译码 code 从未获得译码机会”的问题。

### 1.2 命题 5.2：最终结果等价性

在 5.1 成立的前提下，Buffered FIFO 的每个 code 是否使用了与参考方案一致的：

- EarlyStop 判定；
- channel LLR、先验 / history / extrinsic 输入；
- HISO 或 SISO 动作；
- 写回位置、写回数值与写回顺序；
- 最终 LLR 和信息比特判决。

只有这些条件都成立，才可以说 Buffered FIFO 在最终 BER 意义上等价于“Level 5 和 Level 6 各自具有充足 32 HISO + 32 SISO 资源”的参考方案。

> `所有非 EarlyStop code 都被至少译码一次` 是命题 5.2 的必要条件，但不是充分条件。若一个 code 在两条路径中使用的动作不同，或其 history / extrinsic 输入不同，最终 LLR 和 BER 仍可能不同。

---

## 2. 当前已经确认的事实

固定实验条件：

```text
Eb/N0         = 3.09 dB
bitgen seed   = 20260319
channel seed  = 3182026
R_buf         = 138 block rows
HISO/SISO     = 8 / 8（Buffered FIFO）
```

运行日志：

```text
data/run_20260819-233407_single.log
```

code 生命周期记录：

```text
data/level56_schedule/ofec_single_ebn03.09_rbuf138_level56_schedule_codes.csv
```

对该 CSV 按唯一键 `(batch_id, code)` 聚合后的结果为：

| 检查项 | 结果 |
|---|---:|
| 已进入并实际被服务的 batch | 8,425 |
| 已观察 code 数 | 539,200 = 8,425 × 64 |
| 初始 EarlyStop code | 497,713 |
| 初始非 EarlyStop code | 41,487 |
| 非 EarlyStop 且至少一次 HISO/SISO | 41,487 |
| 非 EarlyStop 但从未 HISO/SISO | **0** |
| 最终记录仍为 pending 的 code | **0** |
| forced-evicted code | **0** |
| `AlreadyDecoded != 1` 的最终 code | **0** |
| `writeback_complete != 1` 的最终 code | **0** |
| EarlyStop 初始/最终分类不一致 | **0** |

因此，对于这 8,425 个已经完成服务的 batch，命题 5.1 已经得到直接证据支持：没有发现任何“初始非 EarlyStop、但整个生命周期从未得到 HISO/SISO”的 code。

### 2.1 当前证据的边界：有限帧尾部没有 drain

本次输入共到达 8,449 个 batch，但只完成服务 8,425 个 batch，差额为 24 个 batch。这些 batch 位于有限帧尾部，当前运行结束时仍在 FIFO 中，未生成完整的 code 生命周期记录。

当前 BER 实现通过裁剪帧尾两个 window 避开该未 drain 尾部。因此：

- 可以确认：已完成、已进入 BER 有效比较区间的服务部分没有遗漏译码机会；
- 不能确认：这一次有限输入中的全部 8,449 个 arrival batch 都完成了 5.1；
- 不能把当前结果称为“整帧全 code 覆盖验证”。

后续验证必须增加显式 tail drain，或明确只对与 BER 有效区间严格对应的 batch 集合做检查。

---

## 3. 参考基线：先使用现有的充足资源 shared 路径

第 2 步的第一参考方案采用现有、已经可运行的配置：

```text
Level5/6 shared              = ON
LEVEL56_SCHEDULE_MODE        = GlobalPriority
shared HISO/SISO             = 32 / 32
temporal lookahead           = OFF
buffered FIFO                = OFF
```

以下称为：

```text
充足资源 shared 参考（shared-32/32 reference）
```

它是 5.2 的合适第一参考，原因是：

- 它复用同一套 Level 5/6 shared 框架、EarlyStopAction、hybrid 分类、HISO/SISO 执行器和 history 写回逻辑；
- 因此与 Buffered FIFO 的主要差别集中在：`8/8` 资源受限后产生 pending、FIFO 延后服务、以及跨时刻补解；
- 若该配置在一个样本中所有非 EarlyStop code 均在到达时获得 HISO/SISO，则它正好表示“没有因为 8/8 burst 而漏解 code”的即时服务参考；
- 现有 `ofec_single` 已验证此路径可直接运行。例如在 `Eb/N0=3.07 dB` 的既有日志中，`shared=ON, 32/32, GlobalPriority` 的 Post-FEC BER 为 0。

### 3.1 它与“Level 5、Level 6 各自独立 32/32”的关系

`shared=ON, 32/32` 是一个 **64 code 共用 32 HISO + 32 SISO** 的资源池；严格意义上，它不等同于：

```text
Level 5 自己有 32 HISO + 32 SISO
Level 6 自己有 32 HISO + 32 SISO
```

二者只会在一个极端情况下不同：例如 64 个 code 中出现超过 32 个只能走 SISO 的 `HardFail` code 时，shared-32/32 的 32 个 SISO 可能仍不足；而两级各自 32 SISO 时，总计可容纳 64 个 SISO-only code。

但这不妨碍 shared-32/32 用作当前 5.2 的第一参考。只要逐 code trace 证明在选定样本中：

```text
shared-32/32 reference 的 non-EarlyStop Unscheduled = 0
```

它就已经具有本次验证需要的性质：参考路径不存在因五六级资源不足而漏解的 code。

`LEVEL56_SHARED_ENABLE=false` 的旧独立 tile 路径仍存在，可在未来作为第二层、更严格的“每级独立资源”参考；但它不是当前 5.2 的前置条件，也不需要为开始第 2 步而先重写。

### 3.2 shared-32/32 参考配置验收条件

建立 A/B 对比前，必须在日志和 code trace 中确认参考路径满足：

- 每个非 EarlyStop code 都有确定的 HISO 或 SISO 最终动作；
- `Unscheduled=0`，不存在资源不足漏解；
- 两级均执行和 Buffered FIFO 相同的 EarlyStop 条件及 `EarlyStopAction`；
- Level 1–4、信道、量化、Chase、alpha/beta、hybrid 参数完全一致；
- 两边使用相同的 BER 比较区间。

---

## 4. 验证总原则

验证顺序必须由内向外：

```text
输入 / 分类一致性
    → code 覆盖与动作一致性
    → 写回一致性
    → 最终 LLR 一致性
    → 最终 bit 判决与 BER 一致性
```

不能从“BER 有差异”直接猜测是 FIFO、EarlyStop、history 还是物理地址映射导致。

当首次出现不一致时，应报告第一个不一致的稳定 code ID，而不是只报告一个总数。

稳定 ID 定义为：

```text
(frame_id, batch_id, level, local_code)
```

其中：

- `batch_id` 为 FIFO arrival 顺序；
- `level ∈ {5, 6}`；
- `local_code ∈ [0, 31]`；
- 同时记录该 code 在 arrival 时的 `source_global_row`，便于回到原始矩阵位置。

---

## 5. 第一阶段：5.1 覆盖性验证

### 5.1 测试目标

对每个被比较的 batch/code，严格确认是否获得了正确且完整的一次处理机会。

### 5.2 每 code 必须记录的不可变字段

在 batch 第一次分类时记录，并在后续 retry 中原样保留：

| 字段 | 含义 |
|---|---|
| `batch_id` | FIFO batch 标识 |
| `level`, `local_code` | Level 5/6 与 0–31 行号 |
| `source_global_row_at_arrival` | arrival 时对应的全局逻辑行 |
| `early_stop_initial` | 初始 EarlyStop 判定 |
| `hybrid_class_initial` | 初始 hybrid 分类 |
| `eligibility_initial` | HISO-only / SISO-only / flexible |
| `channel_llr_hash_at_arrival` | 256 个 channel LLR 的稳定 hash |
| `decoder_input_hash_at_arrival` | 进入 Level 5/6 前 256 个输入 LLR 的稳定 hash |

其中 hash 用于先定位差异。只有 hash 不一致时，再针对该 stable ID 导出完整 256 个 LLR。

### 5.3 每次服务尝试必须记录的字段

| 字段 | 含义 |
|---|---|
| `service_t` | 本次 FIFO 服务时刻 |
| `attempt_index` | 同一 code 的第几次参与调度 |
| `action` | EarlyStopAction / HisoDecode / SisoDecode / Unscheduled |
| `scheduled` | 本轮是否实际获得 HISO/SISO |
| `pending_after_attempt` | 本次后是否仍 pending |
| `already_decoded_after_attempt` | 本次后是否已经完成 |
| `writeback_complete_after_attempt` | 写回是否完成 |
| `output_llr_hash` | 产生结果后的 256 LLR hash；未产生时为空 |
| `history_write_hash` | 如该行写入末级 history，则记录写入值 hash |

### 5.4 5.1 的必须通过断言

对完成 drain 后的每个 batch：

```text
assert(batch.code_count == 64)

for every code:
  if early_stop_initial:
      assert(executed(EarlyStopAction))
      assert(already_decoded && writeback_complete)
  else:
      assert(executed(HisoDecode) || executed(SisoDecode))
      assert(already_decoded && writeback_complete)

assert(final_pending_count == 0)
assert(final_forced_evicted_count == 0)
assert(final_action_failed_count == 0)
```

同时输出以下汇总量：

```text
N_total
N_early_stop
N_non_early_stop
N_non_early_stop_ever_hiso
N_non_early_stop_ever_siso
N_non_early_stop_never_decoded
N_retry_code
N_retry_attempt
N_forced_evicted
N_action_failed
N_writeback_incomplete
```

`N_non_early_stop_never_decoded` 必须为 0，才可说“未获得译码机会问题被消除”。

### 5.5 有限帧 tail drain

为了让 5.1 覆盖所有 arrival batch，验证模式应提供：

```text
输入帧结束
    → 停止注入新 batch
    → 继续执行 FIFO 服务
    → 直到 FIFO 为空
```

这不是改变稳态硬件吞吐的方案声明，而是有限长度软件仿真的收尾操作。它的目的仅是：

- 获取尾部每个 batch 的最终生命周期；
- 让全帧 Post-FEC BER 有明确定义；
- 避免把“未 drain 的尾部”误判为算法差异或缓冲区不足。

验证报告必须把稳态期间的 `forced_evicted` 与 drain 期间的处理结果分开统计。

---

## 6. 第二阶段：5.2 等价性验证

### 6.1 目标

在同一帧、同一信道、同一前四级结果下，比较：

```text
A：Level 5/6 shared-32/32 充足资源参考方案
B：Level 5/6 共享 8 HISO + 8 SISO + Buffered FIFO
```

对 B，只接受通过第 5 节断言、且 `forced_evicted=0` 的运行。

### 6.2 第一层比较：分类输入是否一致

对每一个 stable ID 比较：

| 比较项 | 通过条件 |
|---|---|
| `source_global_row_at_arrival` | 一致 |
| channel LLR hash | 一致 |
| decoder input hash | 一致 |
| EarlyStop 判定 | 一致 |
| hybrid 分类 | 一致 |
| resource eligibility | 一致 |

若这里已经不一致，则问题发生在 Level 5/6 调度器之前，应优先检查：

- 连续 SRAM 逻辑行到全局矩阵行的映射；
- FIFO 等待时 Level 1–4 是否意外修改同一逻辑数据；
- history 输入是否在到达和实际服务之间发生变化；
- batch `tile_top5/tile_top6` 是否仍指向正确逻辑行。

### 6.3 第二层比较：动作是否一致

对于 EarlyStop code，两边必须一致为：

```text
EarlyStopAction
```

对于非 EarlyStop code，需要先记录两边动作对：

```text
(reference_action, buffered_action)
```

并统计：

```text
HISO → HISO
SISO → SISO
HISO → SISO
SISO → HISO
```

后两类并不自动意味着实现错误，但它们是最终 LLR/BER 可能不同的直接候选原因。只有确认 HISO 和 SISO 对这些 code 在当前实现中输出等价，才能忽略这种动作差异。

### 6.4 第三层比较：写回和 history 是否一致

对每个产生结果的 code，记录：

| 字段 | 用途 |
|---|---|
| `writeback_global_row` | 检查写回逻辑地址 |
| `output_llr_hash` | 检查 256 维输出值 |
| `history_input_hash` | 检查进入末级的 history 输入 |
| `history_write_hash` | 检查末级 history 写回 |
| `writeback_order_index` | 检查同一逻辑位置相关的先后关系 |

比较规则：

```text
若 action 相同：output_llr_hash 必须相同。
若 action 不同：记录为 action-induced mismatch，不可直接归咎于地址映射。
若输入 hash 已不同：停止向后归因，先修正输入不一致。
```

### 6.5 第四层比较：最终结果

在完成 tail drain 后，比较：

- 全帧最终 LLR hash；
- 每一信息 bit 的 hard decision；
- 全帧 Post-FEC error positions；
- BER 的错误数和分母。

最有诊断价值的输出是：

```text
first_mismatch_bit
  → global row / col
  → batch_id / level / local_code
  → reference action / buffered action
  → 输入 hash、输出 hash、history hash
```

这样能把“最终 BER 不同”回溯到一个具体 code，而不是只看总错误数猜原因。

---

## 7. 推荐测试顺序

### Test A：单元级覆盖性测试

人工构造一个 batch：

```text
64 code 中有超过一次 8/8 服务能力的非 EarlyStop 请求
```

验证其经历多次 retry 后：

- 每个 non-EarlyStop code 恰好完成一次 HISO/SISO；
- 每个 EarlyStop code 恰好完成一次 EarlyStopAction；
- 已完成 code 不再次进入资源分配；
- 输出和写回次数均为一次。

该测试不依赖 AWGN，也不依赖 BER，可定位 FIFO 状态机错误。

### Test B：真实 trace 的覆盖性回放

选取 `Eb/N0=3.07 dB` 的实际 burst 区间，冻结每个 batch 的初始分类和 LLR 输入，回放 FIFO 调度。

扫描不同 `R_buf`，验证：

```text
R_buf 增大
  → forced_evicted 单调不增（同一冻结 trace）
  → non-EarlyStop never-decoded 单调不增
```

该测试把“调度/缓冲容量问题”与“解码输出改变分类”的反馈分开。

### Test C：端到端 5.1 验证

使用当前真实链路：

```text
Eb/N0=3.09 dB, R_buf=138
```

加入 tail drain，运行第 5 节断言。通过后才说明真实链路中 5.1 对整帧成立。

### Test D：端到端 5.2 A/B 对照

同一帧、同一信道，运行 shared-32/32 充足资源参考方案与通过 Test C 的 Buffered FIFO 方案。按第 6 节四层比较，定位第一个差异。

### Test E：多种子 / 多 Eb/N0 稳健性

在多个 `bitgen_seed`、`channel_seed` 以及 3.07–3.10 dB 范围内重复 Test C/D。

输出的不是单个“最小 R_buf”，而是：

```text
每个 Eb/N0、每个 seed 的 max FIFO depth
每个 Eb/N0 的最坏 required R_buf
forced=0 的概率 / 最坏值
BER 与 32/32 基线的差异分布
```

---

## 8. 已实现的只读观测与开关回归（2026-08-24）

为开始 5.2 的逐 code 对齐，已经加入默认关闭的
`LEVEL56_EQUIVALENCE_OBSERVATION_ENABLE`。`ofec_single` 可通过：

```bash
LEVEL56_EQUIVALENCE_OBSERVATION=1
```

打开。它单独输出：

```text
data/level56_observation/<label>_level56_equivalence_observation.csv
```

每一行是一次实际 Level 5/6 shared 服务中的一个 code。稳定标识为：

```text
(invocation, buffered_t, buffered_batch, code,
 source_level, source_local_row, source_global_row)
```

除既有的 early-stop / hybrid class / HISO-SISO 计划与实际动作外，CSV
还记录以下**只读** FNV-1a-64 指纹：

- `channel_input_hash`：该 code 的 channel 输入；
- `decoder_input_hash`：实际进入 common decoder 的 LIN 输入；
- `decoder_output_hash`：EarlyStopAction / HISO / SISO 后 decoder 输出；
- `tile_row_hash_after_writeback`：该 code 对应 tile 行在 writeback 后的值；
- `level6_history_hash_after_writeback`：仅 Level 6，按原 writeback 映射
  取回的 `last_tile_history_accum` 位置的值。

观测逻辑只在已有计算、分类、调度、执行、writeback 和 history 写回都
完成后读取结果；不会改写 `work_llr`、`last_tile_history_accum`、dispatch、
FIFO 状态或随机数状态。

已完成一次开关回归，固定配置为：

```text
Eb/N0=3.09 dB
R_buf=138
bitgen/channel seed=20260319/3182026
shared FIFO=ON, Group4 8/8
```

对比结果：

| 项目 | observation OFF | observation ON | 结果 |
|---|---:|---:|---|
| Pre-FEC BER | 0.0217743 (622759/28600704) | 0.0217743 (622759/28600704) | 一致 |
| Quantized-hard BER | 0.0218327 (624432/28600704) | 0.0218327 (624432/28600704) | 一致 |
| Post-FEC BER | 2.79713e-07 (8/28600704) | 2.79713e-07 (8/28600704) | 一致 |
| 全部 Post-FEC error positions | 8 个位置 | 同一 8 个位置 | 文件逐字节一致 |
| FIFO 汇总 | completed=8425, forced=0, max depth=70 | 相同 | 一致 |
| 调度 calls / branch | 8554 / [150,8230,145,29] | 相同 | 一致 |

带 history 指纹的 ON 运行导出 `547456 = 8554 × 64` 条观测记录：

```text
decoder_output_available = 539200
level6_history_available = 273728
missing_required_hash = 0
```

这里 `539200 = 8425 × 64`：恰好对应已经退休的 8425 个 batch；帧尾尚在
FIFO 的 24 个 batch 没有完整服务，因而不会产生 decoder-output/history
指纹。`Unscheduled` 的 8256 条尝试记录中，8017 条是已经完成 code 在后续
retry 时被 `AlreadyDecoded` 屏蔽；其余 239 条是该次服务后仍 pending 的 code。
这与此前的生命周期聚合统计一致，不表示有 code 丢失。

因此，观测开关本身已经以端到端固定样本证明不改变解码结果。下一轮 5.2
对照可直接拿这份 CSV 与 shared-32/32 参考 CSV 进行逐 code 对齐。

## 9. 5.2 第一层 B/C/D 控制实验（2026-08-24）

### 9.1 三个控制组

在与第 8 节相同的 bit/channel seed、Eb/N0、量化、early-stop action 和
hybrid 分类参数下，完成以下三组端到端运行：

| 组别 | shared | FIFO | HISO/SISO | scheduler | R_buf | 用途 |
|---|---|---|---:|---|---:|---|
| B | ON | OFF | 32/32 | GlobalPriority | N/A | 即时充足资源参考 |
| D | ON | ON | 32/32 | GlobalPriority | 138 | FIFO 读/写/延后语义控制组 |
| C | ON | ON | 8/8 | Group4 multiround | 138 | 新方案实际配置 |

D 是必要控制组：它将 FIFO 存储语义与 Group4 的 HISO/SISO 动作选择分开。
Group4 保持原有固定 8-entry、8/8 约束；GlobalPriority 则由原 common
scheduler 按配置的 32/32 容量调度。FIFO 状态机、physical-SRAM read/write、
固定吞吐和 forced-eviction 逻辑均未改变。

### 9.2 端到端结果

| 组别 | Post-FEC BER | error positions | FIFO 完成 / forced / max depth |
|---|---:|---|---|
| B | 0 / 28600704 | 空 | 不适用 |
| D | 0 / 28600704 | 空，和 B 同一空文件 hash | 8449 / 0 / 1 |
| C | 8 / 28600704 = 2.79713e-07 | 8 个位置 | 8425 / 0 / 70 |

B 与 D 的最终 hard decision 完全一致；C 的 8 个错误在 clean rebuild 后再次
复现，和此前 C 运行的 error-position 文件逐字节一致。

### 9.3 B 与 D 的逐 code 对齐结论

B/D 的观测 CSV 各有 `8449 × 64 = 540736` 条记录。按 `(batch_id, code)`
对齐，以下字段都是 **0 个差异**：source global row、early-stop hit、
hybrid class / eligibility、final action、HISO/SISO plan 与 core assignment、
channel_input_hash、decoder_input_hash、decoder_output_hash、
tile_row_hash_after_writeback、level6_history_hash_after_writeback。

唯一不同是预期的 FIFO 内部状态字段：D 的 code 服务后标为
`AlreadyDecoded/writeback_complete`，B 不使用该状态；它不对应任何 LLR、
history 或 BER 差异。

因此，在这个固定样本下可以确认：FIFO 的连续物理 SRAM 读取/写回、FIFO
退休和 history 写回语义本身，不是 C 相比 B 出现 8 个 Post-FEC 错误的原因。

### 9.4 B 与 C：第一个真实分叉

B/C 首次按 `(batch_id, code)` 对齐时，第一个差异为：

```text
batch B50, Level 5, local code 4, global row 2275
channel_input_hash : 相同
decoder_input_hash : 相同
B final_action     : SisoDecode
C final_action     : HisoDecode
```

这是同一个输入下的动作选择不同，发生在任何后续 LLR/history 分叉之前。
首个服务尝试的动作差异总计为：

```text
B SisoDecode -> C HisoDecode       7887
B SisoDecode -> C Unscheduled       237
B SisoDecode -> C EarlyStopAction     34
B EarlyStopAction -> C HisoDecode      9
B EarlyStopAction -> C SisoDecode      12
B EarlyStopAction -> C Unscheduled      2
```

随后 B/C 的 `decoder_input_hash` 出现大面积差异，这是首个 action 差异写回后
的因果后果，不应误判为 FIFO 延后读到了错误输入。

动作总量也说明三组不是同一调度策略：

```text
B/D: EarlyStopAction=499108, SisoDecode=41628, HisoDecode=0, Unscheduled=0
C  : EarlyStopAction=497713, SisoDecode=33575, HisoDecode=7912,
     Unscheduled attempt=8256（其中 8017 为 AlreadyDecoded 屏蔽，239 为 pending）
```

所以当前 C 对 B 的 BER 差异，不能解释成“所有 code 获得机会以后，FIFO 延后仍然
损失性能”。C 与 B 同时改变了资源容量和调度策略：Group4 8/8 把一部分 B 中走
SISO 的 flexible code 改成 HISO，并在部分 code 上延后。下一层验证需要保持
Group4 action policy 的容量控制组，才能继续把“动作类型差异”与“延后补解差异”
分开。

## 10. Group4 SISO-only 控制实验：0/8 FIFO 与 0/32 参考（2026-08-24）

### 10.1 目的与两组配置

第 9 节说明原 `8 HISO / 8 SISO` Group4 与 32/32 参考之间首先存在
`SisoDecode → HisoDecode` 的动作差异。为把这个因素完全移除，本节新增了一个
**显式 SISO-only Group4** 模式：

```text
LEVEL56_SHARED_HISO_ACTIVE = 0
LEVEL56_SHARED_SISO_ACTIVE = 8
LEVEL56_SCHEDULE_MODE      = Group4LoadSortedMultiround
```

它不是“把原 8/8 的 HISO 参数改成零”。其定义为：每个获得 Group4 entry 的组最多
选一个 code 执行 SISO，绝不产生 HISO；同组其余 non-EarlyStop code 留作 pending，
在后续 FIFO 时刻重新进入原 Group4 调度器。

在固定样本上比较：

```text
Eb/N0              = 3.09 dB
bit/channel seed   = 20260319 / 3182026
early-stop、hybrid、量化、decoder 均保持相同
```

| 组别 | FIFO | HISO/SISO | scheduler | R_buf | 帧尾 drain | 用途 |
|---|---|---:|---|---:|---|---|
| E（参考） | OFF | 0/32 | GlobalPriority | N/A | N/A | 每个 non-EarlyStop code 立即获得 SISO |
| K（实验） | ON | 0/8 | Group4 SISO-only | 1600 | ON | 仅保留 FIFO 延后与 Group4 entry 限制 |

`R_buf=1600` 的正常到达阶段从未顶到边界；这点由
`forced_evicted=0` 和 `min_S_t=124` 共同确认。

### 10.2 为什么必须加帧尾 drain

正常固定吞吐仿真在最后一个输入 batch 到达后即停止。对本实验的 `R_buf=1600`，此时：

```text
arrival batch = 8449
completed     = 7748
FIFO tail     = 701 batch
forced        = 0
```

这里的 701 个 batch 是“尚未完成”，不是“被 forced eviction”。旧 BER 的固定尾部只
裁两个 window，远不足以排除这段未完成 FIFO tail；因此不能用它与立即完成的参考组
讨论最终 BER。

为此加入默认关闭的验证开关：

```bash
LEVEL56_BUFFERED_FIFO_DRAIN_AT_FRAME_END=1
```

其语义严格限定为：所有正常 arrival time 已结束后，**不引入新 batch、不执行
Level 1--4**，只按原 FIFO 规则继续服务队首，直到 FIFO 清空。它不改变任何帧内时刻的
调度顺序，也不是在线固定吞吐方案的一部分。

### 10.3 第一次带 drain 的观察：定位 history 覆盖

在 SISO-only、`R_buf=1600`、drain ON 下，第一次运行已满足：

```text
completed=8449, forced_evicted=0
HISO=0, SISO=41628
```

但得到 `58 / 28600704` 个 Post-FEC errors，而 E 为 0。逐 code 对齐结果是：

```text
41628 个 SISO：decoder_input_hash、decoder_output_hash、
              tile_row_hash_after_writeback 全部 0 差异
```

第一个差异出现在 Level 6 `last_tile_history_accum`。根因是旧
`writeback_tile()` 对“本次没有 produced”的行仍做 history prior 透传；FIFO batch
多次服务时，pending/AlreadyDecoded 行会再次经过该分支，从而覆盖同一物理位置上此前
已经写好的最终 Level 6 history。参考组每个 batch 一次完成，不会产生这个中间覆盖。

这不符合新方案已确定的语义：**没有新结果，不写回**。修正为：

```text
普通（非 FIFO）路径：保持原 history prior 透传行为；
Buffered FIFO 路径：仅本次 produced 的 code 更新 Level 6 history；
                    pending / AlreadyDecoded 不覆盖既有 history。
```

### 10.4 正式结果：本固定样本达到与 0/32 参考相同的 BER

在上述 history 修正后，重新 clean build、通过既有 shared scheduler regression check，
并重跑完全相同的 K 配置。结果如下：

| 指标 | E：FIFO OFF, 0/32 | K：FIFO ON, 0/8 Group4, R_buf=1600, drain ON |
|---|---:|---:|
| Post-FEC BER | `0 / 28600704` | `0 / 28600704` |
| Post-FEC error positions | 空 | 空，文件逐字节一致 |
| HISO 调用 | 0 | 0 |
| SISO 调用 | 41628 | 41628 |
| completed batch | 8449 | 8449 |
| forced evicted batch | N/A | 0 |
| full-EarlyStop batch | N/A | 150 |
| 最大 FIFO 深度 | N/A | 702 |
| `min_S_t` | N/A | 124 |

对 K 的 8449 × 64 个 code 生命周期聚合：

```text
initial non-EarlyStop = 41628
获得至少一次 SISO       = 41628
从未获得 SISO           = 0
HISO                    = 0
initial EarlyStop       = 499108
获得 EarlyStopAction    = 499108
forced eviction         = 0
```

因此，在**这个固定 Eb/N0/seed 样本**中，已经实测支持下述结论：

> 若缓冲区足够长、无 forced eviction、帧尾 FIFO 完整 drain，且 pending /
> AlreadyDecoded code 不以“无新结果”的形式覆盖 Level 6 history，那么 `0 HISO +
> 8 SISO` 的 FIFO Group4 可以让全部需要软译码的 code 最终各获得一次 SISO，并达到与
> `FIFO OFF + 0 HISO + 32 SISO` 参考组完全相同的最终 hard decision / BER。

这不是对所有 Eb/N0、所有 seed 的数学证明；它是一个严格控制的端到端样本验证。下一步应
在多个 seed 和 3.07--3.10 dB 点重复本节 E/K 对照，并记录达到 `forced=0` 所需的最坏
`R_buf` 和帧尾 drain 长度。
