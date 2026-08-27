# Level 5/6 Buffered FIFO：不同解码配置下 batch 服务次数汇总

> 生成时间：2026-08-26。服务次数定义为 `ordinary_batch` 在时间 CSV 中出现的次数；`FullEarlyStop` 不占普通 HISO/SISO 预算，因此普通服务次数为 0。每张表只对已经退休的 batch 统计生命周期分布；未 drain 运行中仍留在 FIFO 的 batch 计入 `censored`，不计入分位数。百分位采用线性插值。

## 字段说明

- `max_pending`：该 batch 被观察到的最大 `pending_before`；比退休后的 `final pending` 更有信息，因为退休时 pending 会被清零。
- `P95 max_pending`：先对每个 batch 取其生命周期内的最大 pending，再对这些 batch 级最大值取 95 分位；它不是某一时刻所有 pending 值的 95 分位。
- `服务4次以上` 是服务次数 `>=4` 的 batch 数。

## 各配置 Pre-FEC / Post-FEC 错误汇总

| 配置 | run_id | Pre-FEC错误数 | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | BER比较总比特数 |
|---|---|---:|---:|---:|---:|---:|
| Eb/N0=3.065 dB; FIFO=ON; R_buf=16000; HISO/SISO=0/8; Group4; entries=default; drain=ON | `20260825_002320` | 1292037 | 0.0220454 | 11 | 1.87688e-07 | 58608000 |
| Eb/N0=3.065 dB; FIFO=ON; R_buf=16000; HISO/SISO=8/8; Group4; entries=default; drain=ON | `20260825_150958` | 1292037 | 0.0220454 | 159 | 2.71294e-06 | 58608000 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=ON | `20260824_212534` | 622759 | 0.0217743 | 0 | 0 | 28600704 |
| Eb/N0=3.1 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_214550` | 619300 | 0.0216533 | 5 | 1.74821e-07 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `20260824_203748` | 622759 | 0.0217743 | 341 | 1.19228e-05 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=512; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `20260824_204440` | 622759 | 0.0217743 | 14649 | 0.00051219 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=1200; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `20260824_205202` | 622759 | 0.0217743 | 41063 | 0.00143573 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `20260824_210302` | 622759 | 0.0217743 | 48920 | 0.00171045 | 28600704 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=ON | `20260824_231958` | 633274 | 0.0221419 | 27449 | 0.000959732 | 28600704 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=3200; HISO/SISO=0/8; Group4; entries=default; drain=ON | `20260824_232620` | 633274 | 0.0221419 | 14834 | 0.000518659 | 28600704 |
| Eb/N0=3.065 dB; FIFO=ON; R_buf=6400; HISO/SISO=0/8; Group4; entries=default; drain=ON | `20260824_235918` | 1292037 | 0.0220454 | 13408 | 0.000228774 | 58608000 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_161506` | 591774 | 0.0221429 | 3431 | 0.00012838 | 26725248 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_155111` | 638406 | 0.0221399 | 3525 | 0.000122247 | 28835136 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_155428` | 638406 | 0.0221399 | 3724 | 0.000129148 | 28835136 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_162031` | 633274 | 0.0221419 | 3487 | 0.00012192 | 28600704 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_163014` | 629688 | 0.0220165 | 544 | 1.90205e-05 | 28600704 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_163337` | 629688 | 0.0220165 | 543 | 1.89855e-05 | 28600704 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_163656` | 629688 | 0.0220165 | 543 | 1.89855e-05 | 28600704 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_164023` | 629688 | 0.0220165 | 544 | 1.90205e-05 | 28600704 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_224406` | 629688 | 0.0220165 | 544 | 1.90205e-05 | 28600704 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=80; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_230543` | 629688 | 0.0220165 | 520 | 1.81814e-05 | 28600704 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=96; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_230948` | 629688 | 0.0220165 | 503 | 1.7587e-05 | 28600704 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_164348` | 626242 | 0.021896 | 153 | 5.34952e-06 | 28600704 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_164723` | 626242 | 0.021896 | 126 | 4.40549e-06 | 28600704 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_165038` | 626242 | 0.021896 | 128 | 4.47541e-06 | 28600704 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_165349` | 626242 | 0.021896 | 128 | 4.47541e-06 | 28600704 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_231259` | 626242 | 0.021896 | 128 | 4.47541e-06 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_231647` | 622759 | 0.0217743 | 34 | 1.18878e-06 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=80; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_232033` | 622759 | 0.0217743 | 30 | 1.04893e-06 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=128; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_232417` | 622759 | 0.0217743 | 9 | 3.14678e-07 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=136; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_233726` | 622759 | 0.0217743 | 8 | 2.79713e-07 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260824_172806` | 622759 | 0.0217743 | 8 | 2.79713e-07 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=140; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_233054` | 622759 | 0.0217743 | 8 | 2.79713e-07 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=144; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_232736` | 622759 | 0.0217743 | 8 | 2.79713e-07 | 28600704 |
| Eb/N0=3.1 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `20260819_200826` | 619300 | 0.0216533 | 6 | 2.09785e-07 | 28600704 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=32/32; GlobalPriority; entries=default; drain=OFF | `20260824_194542` | 622759 | 0.0217743 | 0 | 0 | 28600704 |

## BER 随 R_buf / HISO-SISO 的变化分析

> 下面的趋势只使用上方 BER 汇总表中实际存在的运行。不同 Eb/N0、帧长、drain 状态或随机输入不混合计算；`R_buf` 变化趋势优先在固定 Eb/N0、资源和调度模式下观察。

### 1. 8/8 配置下，Post-FEC BER 随 R_buf 的变化

这里的 8/8 指 `HISO/SISO=8/8`，调度模式为 Group4。

#### Eb/N0=3.06 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.06 | 0 | OFF | 0.0221429 | 3431 | 0.00012838 | 26725248 |
| 3.06 | 8 | OFF | 0.0221399 | 3525 | 0.000122247 | 28835136 |
| 3.06 | 16 | OFF | 0.0221399 | 3724 | 0.000129148 | 28835136 |
| 3.06 | 32 | OFF | 0.0221419 | 3487 | 0.00012192 | 28600704 |

结论：Post-FEC BER 随 R_buf 不是严格单调，存在运行噪声或其他状态因素。
> 注意：该组包含 drain=OFF 运行；BER 可能受到尾部未完成/ForcedEvicted 影响，不能只归因于 R_buf。

#### Eb/N0=3.065 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.065 | 16000 | ON | 0.0220454 | 159 | 2.71294e-06 | 58608000 |

#### Eb/N0=3.07 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.07 | 0 | OFF | 0.0220165 | 544 | 1.90205e-05 | 28600704 |
| 3.07 | 8 | OFF | 0.0220165 | 543 | 1.89855e-05 | 28600704 |
| 3.07 | 16 | OFF | 0.0220165 | 543 | 1.89855e-05 | 28600704 |
| 3.07 | 32 | OFF | 0.0220165 | 544 | 1.90205e-05 | 28600704 |
| 3.07 | 64 | OFF | 0.0220165 | 544 | 1.90205e-05 | 28600704 |
| 3.07 | 80 | OFF | 0.0220165 | 520 | 1.81814e-05 | 28600704 |
| 3.07 | 96 | OFF | 0.0220165 | 503 | 1.7587e-05 | 28600704 |

结论：Post-FEC BER 随 R_buf 不是严格单调，存在运行噪声或其他状态因素。
> 注意：该组包含 drain=OFF 运行；BER 可能受到尾部未完成/ForcedEvicted 影响，不能只归因于 R_buf。

#### Eb/N0=3.08 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.08 | 0 | OFF | 0.021896 | 153 | 5.34952e-06 | 28600704 |
| 3.08 | 8 | OFF | 0.021896 | 126 | 4.40549e-06 | 28600704 |
| 3.08 | 16 | OFF | 0.021896 | 128 | 4.47541e-06 | 28600704 |
| 3.08 | 32 | OFF | 0.021896 | 128 | 4.47541e-06 | 28600704 |
| 3.08 | 64 | OFF | 0.021896 | 128 | 4.47541e-06 | 28600704 |

结论：Post-FEC BER 随 R_buf 不是严格单调，存在运行噪声或其他状态因素。
> 注意：该组包含 drain=OFF 运行；BER 可能受到尾部未完成/ForcedEvicted 影响，不能只归因于 R_buf。

#### Eb/N0=3.09 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.09 | 64 | OFF | 0.0217743 | 34 | 1.18878e-06 | 28600704 |
| 3.09 | 80 | OFF | 0.0217743 | 30 | 1.04893e-06 | 28600704 |
| 3.09 | 128 | OFF | 0.0217743 | 9 | 3.14678e-07 | 28600704 |
| 3.09 | 136 | OFF | 0.0217743 | 8 | 2.79713e-07 | 28600704 |
| 3.09 | 138 | OFF | 0.0217743 | 8 | 2.79713e-07 | 28600704 |
| 3.09 | 140 | OFF | 0.0217743 | 8 | 2.79713e-07 | 28600704 |
| 3.09 | 144 | OFF | 0.0217743 | 8 | 2.79713e-07 | 28600704 |

结论：Post-FEC BER 随 R_buf 增大整体不升（单调下降或持平）。
> 注意：该组包含 drain=OFF 运行；BER 可能受到尾部未完成/ForcedEvicted 影响，不能只归因于 R_buf。

#### Eb/N0=3.1 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.1 | 32 | OFF | 0.0216533 | 6 | 2.09785e-07 | 28600704 |
| 3.1 | 64 | OFF | 0.0216533 | 5 | 1.74821e-07 | 28600704 |

结论：Post-FEC BER 随 R_buf 增大整体不升（单调下降或持平）。
> 注意：该组包含 drain=OFF 运行；BER 可能受到尾部未完成/ForcedEvicted 影响，不能只归因于 R_buf。

### 2. 0/8 配置下，Post-FEC BER 随 R_buf 的变化

这里的 0/8 指 `HISO/SISO=0/8`，即纯 SISO 的受限 FIFO 配置。

#### Eb/N0=3.06 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.06 | 1600 | ON | 0.0221419 | 27449 | 0.000959732 | 28600704 |
| 3.06 | 3200 | ON | 0.0221419 | 14834 | 0.000518659 | 28600704 |

结论：Post-FEC BER 随 R_buf 增大整体不升（单调下降或持平）。

#### Eb/N0=3.065 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.065 | 6400 | ON | 0.0220454 | 13408 | 0.000228774 | 58608000 |
| 3.065 | 16000 | ON | 0.0220454 | 11 | 1.87688e-07 | 58608000 |

结论：Post-FEC BER 随 R_buf 增大整体不升（单调下降或持平）。

#### Eb/N0=3.09 dB

| Eb/N0 (dB) | R_buf | drain | Pre-FEC BER | Post-FEC错误数 | Post-FEC BER | 比较总比特数 |
|---:|---:|---|---:|---:|---:|---:|
| 3.09 | 138 | OFF | 0.0217743 | 341 | 1.19228e-05 | 28600704 |
| 3.09 | 512 | OFF | 0.0217743 | 14649 | 0.00051219 | 28600704 |
| 3.09 | 1200 | OFF | 0.0217743 | 41063 | 0.00143573 | 28600704 |
| 3.09 | 1600 | ON | 0.0217743 | 0 | 0 | 28600704 |
| 3.09 | 1600 | OFF | 0.0217743 | 48920 | 0.00171045 | 28600704 |

结论：Post-FEC BER 随 R_buf 不是严格单调，存在明显非单调变化。
> 注意：该组包含 drain=OFF 运行；请结合 censored/ForcedEvicted 数量阅读。

### 3. 同 R_buf 下，8/8 与 0/8 的 Post-FEC BER 对比

只列出当前数据中 Eb/N0、帧长/运行样本和 R_buf 都能对齐的配对；没有对应数据的 R_buf 不做推断。

| Eb/N0 (dB) | R_buf | 8/8 Pre-FEC BER | 0/8 Pre-FEC BER | 8/8 Post-FEC BER | 0/8 Post-FEC BER | 8/8错误数 | 0/8错误数 | 说明 |
|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 3.065 | 16000 | 0.0220454 | 0.0220454 | 2.71294e-06 | 1.87688e-07 | 159 | 11 | 同 Eb/N0、R_buf、帧长；BER比值(8/8 ÷ 0/8)=14.5 |
| 3.09 | 138 | 0.0217743 | 0.0217743 | 2.79713e-07 | 1.19228e-05 | 8 | 341 | 同 Eb/N0、R_buf、帧长；BER比值(8/8 ÷ 0/8)=0.0235 |

#### 解读

- `R_buf` 增大主要改变 pending 的保存深度和 ForcedEvicted 风险；它不直接改变单个 code 的 HISO/SISO 解码算法。
- 因此在同一资源配置下，BER 随 `R_buf` 下降通常表示 buffer 减少了边界淘汰或尾部未完成影响；若不单调，应结合 drain、ForcedEvicted 和 censored 数量判断。
- 8/8 与 0/8 的差异不能只解释为服务次数差异：HISO 会改变 code 的输出和后续 SRAM/history 输入，因此即使两组都没有 ForcedEvicted，BER 也可能不同。
- 当前同 R_buf 配对数量有限，以上是已有样本的直接比较，不代表对所有随机帧和 Eb/N0 的统计规律。

#### 为什么长帧和短帧的结论方向相反

- 长帧 `Eb/N0=3.065, R_buf=16000, drain=ON` 配对中，两组都是 16897/16897 个 arrival batch 已退休，`ForcedEvicted=0`，因此比较的是完整生命周期。0/8 的 Post-FEC BER 为 `1.87688e-7`，8/8 为 `2.71294e-6`。
- 短帧 `Eb/N0=3.09, R_buf=138, drain=OFF` 配对中，8/8 仍有 24 个 censored batch，0/8 有 67 个 censored batch，并且 0/8 有 855 个 `ForcedEvicted`，8/8 为 0。0/8 的 `1.19228e-5` 与 8/8 的 `2.79713e-7` 不能视为只由 HISO/SISO 导致的差异。
- 短帧的两个运行虽然 Pre-FEC BER 相同，但 buffer 压力和退休路径不同：0/8 没有 HISO 服务，SISO 压力更集中，更容易出现 FIFO 边界淘汰；8/8 的 HISO 服务会改变历史状态写回和后续窗口输入。两种机制叠加后，方向完全可能与长帧不同。
- 另外，长帧和短帧的 `R_buf`、drain 状态、观测长度不同；短帧 Post-FEC 错误数只有 8（8/8）或 341（0/8），低错误数的相对统计不确定性也更大。因此这两个配对应作为两个不同实验解读，不能拼成一条“8/8 一定优于/劣于 0/8”的规律。
- 要判断 HISO/SISO 的稳定影响，建议固定 `Eb/N0`、帧长、R_buf、drain、seeds 和 equivalence-observation 设置，至少各重复多次，并优先使用 `drain=ON`、`ForcedEvicted=0`、`censored=0` 的运行。

> Pre-FEC/Post-FEC 错误数和比较总比特数直接取对应 run log 的 `[RESULT]` 行；`?` 表示该日志没有可解析的 BER 结果。

## 数据源清单

| 配置 | CSV | run_id | arrival | retired | censored | rows |
|---|---|---|---:|---:|---:|---:|
| Eb/N0=3.065 dB; FIFO=ON; R_buf=16000; HISO/SISO=0/8; Group4; entries=default; drain=ON | `debug_L6_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso0_siso8_schedgroup4_level56_buffered_times.csv` | `20260825_002320` | 16897 | 16897 | 0 | 22971 |
| Eb/N0=3.065 dB; FIFO=ON; R_buf=16000; HISO/SISO=8/8; Group4; entries=default; drain=ON | `debug_L6_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso8_siso8_schedgroup4_level56_buffered_times.csv` | `20260825_150958` | 16897 | 16897 | 0 | 19109 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=ON | `debug_L6_ebn03.09_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv` | `20260824_212534` | 8449 | 8449 | 0 | 9184 |
| Eb/N0=3.1 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.1_rbuf64_level56_buffered_times.csv` | `20260819_214550` | 8449 | 8449 | 0 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_fifo1_rbuf138_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv` | `20260824_203748` | 8449 | 8382 | 67 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=512; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_fifo1_rbuf512_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv` | `20260824_204440` | 8449 | 8193 | 256 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=1200; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_fifo1_rbuf1200_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv` | `20260824_205202` | 8449 | 7849 | 600 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_fifo1_rbuf1600_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv` | `20260824_210302` | 8449 | 7748 | 701 | 8449 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=ON | `debug_L6_ebn03.06_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv` | `20260824_231958` | 8449 | 8449 | 0 | 9564 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=3200; HISO/SISO=0/8; Group4; entries=default; drain=ON | `debug_L6_ebn03.06_fifo1_rbuf3200_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv` | `20260824_232620` | 8449 | 8449 | 0 | 10572 |
| Eb/N0=3.065 dB; FIFO=ON; R_buf=6400; HISO/SISO=0/8; Group4; entries=default; drain=ON | `debug_L6_ebn03.065_nbits60014592_fifo1_rbuf6400_drain1_hiso0_siso8_schedgroup4_level56_buffered_times.csv` | `20260824_235918` | 16897 | 16897 | 0 | 21385 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_fifo1_rbuf0_level56_buffered_times.csv` | `20260819_161506` | 8449 | 8449 | 0 | 8449 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_fifo1_rbuf8_level56_buffered_times.csv` | `20260819_155111` | 8449 | 8446 | 3 | 8449 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_fifo1_rbuf16_level56_buffered_times.csv` | `20260819_155428` | 8449 | 8442 | 7 | 8449 |
| Eb/N0=3.06 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_fifo1_rbuf32_level56_buffered_times.csv` | `20260819_162031` | 8449 | 8434 | 15 | 8449 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.07_fifo1_rbuf0_level56_buffered_times.csv` | `20260819_163014` | 8449 | 8449 | 0 | 8449 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.07_fifo1_rbuf8_level56_buffered_times.csv` | `20260819_163337` | 8449 | 8445 | 4 | 8449 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.07_fifo1_rbuf16_level56_buffered_times.csv` | `20260819_163656` | 8449 | 8441 | 8 | 8449 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.07_fifo1_rbuf32_level56_buffered_times.csv` | `20260819_164023` | 8449 | 8433 | 16 | 8449 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.07_rbuf64_level56_buffered_times.csv` | `20260819_224406` | 8449 | 8417 | 32 | 8449 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=80; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.07_rbuf80_level56_buffered_times.csv` | `20260819_230543` | 8449 | 8409 | 40 | 8449 |
| Eb/N0=3.07 dB; FIFO=ON; R_buf=96; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.07_rbuf96_level56_buffered_times.csv` | `20260819_230948` | 8449 | 8401 | 48 | 8449 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.08_fifo1_rbuf0_level56_buffered_times.csv` | `20260819_164348` | 8449 | 8449 | 0 | 8449 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.08_fifo1_rbuf8_level56_buffered_times.csv` | `20260819_164723` | 8449 | 8445 | 4 | 8449 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.08_fifo1_rbuf16_level56_buffered_times.csv` | `20260819_165038` | 8449 | 8441 | 8 | 8449 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.08_fifo1_rbuf32_level56_buffered_times.csv` | `20260819_165349` | 8449 | 8433 | 16 | 8449 |
| Eb/N0=3.08 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.08_rbuf64_level56_buffered_times.csv` | `20260819_231259` | 8449 | 8417 | 32 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_rbuf64_level56_buffered_times.csv` | `20260819_231647` | 8449 | 8449 | 0 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=80; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_rbuf80_level56_buffered_times.csv` | `20260819_232033` | 8449 | 8449 | 0 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=128; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_rbuf128_level56_buffered_times.csv` | `20260819_232417` | 8449 | 8430 | 19 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=136; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_rbuf136_level56_buffered_times.csv` | `20260819_233726` | 8449 | 8426 | 23 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_rbuf138_level56_buffered_times.csv` | `20260824_172806` | 8449 | 8425 | 24 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=140; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_rbuf140_level56_buffered_times.csv` | `20260819_233054` | 8449 | 8425 | 24 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=144; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.09_rbuf144_level56_buffered_times.csv` | `20260819_232736` | 8449 | 8425 | 24 | 8449 |
| Eb/N0=3.1 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF | `debug_L6_ebn03.1_rbuf32_level56_buffered_times.csv` | `20260819_200826` | 8449 | 8449 | 0 | 8449 |
| Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=32/32; GlobalPriority; entries=default; drain=OFF | `debug_L6_ebn03.09_fifo1_rbuf138_hiso32_siso32_schedglobal_eqobs1_level56_buffered_times.csv` | `20260824_194542` | 8449 | 8449 | 0 | 8449 |

## 1. Eb/N0=3.065 dB; FIFO=ON; R_buf=16000; HISO/SISO=0/8; Group4; entries=default; drain=ON

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso0_siso8_schedgroup4_level56_buffered_times.csv`；run_id：`20260825_002320`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 16781 | 5 | 1.36 | 11481 | 4561 | 586 | 131 | 1.45 | 8 |
| FullEarlyStop | 116 | 0 | 0 | 0 | 0 | 0 | 0 | 0.02 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 2. Eb/N0=3.065 dB; FIFO=ON; R_buf=16000; HISO/SISO=8/8; Group4; entries=default; drain=ON

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.065_nbits60014592_fifo1_rbuf16000_drain1_hiso8_siso8_schedgroup4_level56_buffered_times.csv`；run_id：`20260825_150958`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 16781 | 4 | 1.14 | 14668 | 1939 | 155 | 11 | 0.47 | 3 |
| FullEarlyStop | 116 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 3. Eb/N0=3.09 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=ON

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv`；run_id：`20260824_212534`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8299 | 3 | 1.10 | 7466 | 800 | 20 | 0 | 0.25 | 2 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0.04 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 4. Eb/N0=3.1 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.1_rbuf64_level56_buffered_times.csv`；run_id：`20260819_214550`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8260 | 2 | 1.01 | 8203 | 57 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 189 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 5. Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=0/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_fifo1_rbuf138_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv`；run_id：`20260824_203748`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7383 | 2 | 1.02 | 7217 | 166 | 0 | 0 | 0.05 | 0 |
| FullEarlyStop | 144 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 855 | 1 | 1 | 855 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 67 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8382 个已退休 batch。

## 6. Eb/N0=3.09 dB; FIFO=ON; R_buf=512; HISO/SISO=0/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_fifo1_rbuf512_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv`；run_id：`20260824_204440`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7460 | 3 | 1.05 | 7113 | 343 | 4 | 0 | 0.12 | 0 |
| FullEarlyStop | 140 | 0 | 0 | 0 | 0 | 0 | 0 | 0.03 | 0 |
| ForcedEvicted | 593 | 1 | 1 | 593 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 256 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8193 个已退休 batch。

## 7. Eb/N0=3.09 dB; FIFO=ON; R_buf=1200; HISO/SISO=0/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_fifo1_rbuf1200_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv`；run_id：`20260824_205202`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7583 | 3 | 1.09 | 6911 | 652 | 20 | 0 | 0.23 | 2 |
| FullEarlyStop | 137 | 0 | 0 | 0 | 0 | 0 | 0 | 0.04 | 0 |
| ForcedEvicted | 129 | 1 | 1 | 129 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 600 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 7849 个已退休 batch。

## 8. Eb/N0=3.09 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_fifo1_rbuf1600_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv`；run_id：`20260824_210302`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7612 | 3 | 1.10 | 6840 | 752 | 20 | 0 | 0.26 | 2 |
| FullEarlyStop | 136 | 0 | 0 | 0 | 0 | 0 | 0 | 0.04 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 注意：该运行有 701 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 7748 个已退休 batch。

## 9. Eb/N0=3.06 dB; FIFO=ON; R_buf=1600; HISO/SISO=0/8; Group4; entries=default; drain=ON

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.06_fifo1_rbuf1600_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv`；run_id：`20260824_231958`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 5526 | 5 | 1.21 | 4525 | 875 | 112 | 12 | 0.81 | 6 |
| FullEarlyStop | 67 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 2856 | 1 | 1 | 2856 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 10. Eb/N0=3.06 dB; FIFO=ON; R_buf=3200; HISO/SISO=0/8; Group4; entries=default; drain=ON

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.06_fifo1_rbuf3200_drain1_hiso0_siso8_schedgroup4_eqobs1_level56_buffered_times.csv`；run_id：`20260824_232620`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 6543 | 5 | 1.33 | 4755 | 1499 | 230 | 54 | 1.37 | 8 |
| FullEarlyStop | 68 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 1838 | 3 | 1.00 | 1837 | 0 | 1 | 0 | 0.01 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 11. Eb/N0=3.065 dB; FIFO=ON; R_buf=6400; HISO/SISO=0/8; Group4; entries=default; drain=ON

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.065_nbits60014592_fifo1_rbuf6400_drain1_hiso0_siso8_schedgroup4_level56_buffered_times.csv`；run_id：`20260824_235918`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 14889 | 5 | 1.31 | 11005 | 3311 | 453 | 107 | 1.24 | 7 |
| FullEarlyStop | 114 | 0 | 0 | 0 | 0 | 0 | 0 | 0.02 | 0 |
| ForcedEvicted | 1894 | 1 | 1 | 1894 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 12. Eb/N0=3.06 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_fifo1_rbuf0_level56_buffered_times.csv`；run_id：`20260819_161506`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 6737 | 1 | 1 | 6737 | 0 | 0 | 0 | 0 | 0 |
| FullEarlyStop | 69 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 1643 | 1 | 1 | 1643 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 13. Eb/N0=3.06 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_fifo1_rbuf8_level56_buffered_times.csv`；run_id：`20260819_155111`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 6763 | 2 | 1.00 | 6737 | 26 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 69 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 1614 | 1 | 1 | 1614 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 3 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8446 个已退休 batch。

## 14. Eb/N0=3.06 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_fifo1_rbuf16_level56_buffered_times.csv`；run_id：`20260819_155428`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 6764 | 2 | 1.00 | 6733 | 31 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 69 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 1609 | 1 | 1 | 1609 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 7 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8442 个已退休 batch。

## 15. Eb/N0=3.06 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_fifo1_rbuf32_level56_buffered_times.csv`；run_id：`20260819_162031`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 6764 | 2 | 1.01 | 6725 | 39 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 69 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 1601 | 1 | 1 | 1601 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 15 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8434 个已退休 batch。

## 16. Eb/N0=3.07 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.07_fifo1_rbuf0_level56_buffered_times.csv`；run_id：`20260819_163014`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7572 | 1 | 1 | 7572 | 0 | 0 | 0 | 0 | 0 |
| FullEarlyStop | 93 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 784 | 1 | 1 | 784 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 17. Eb/N0=3.07 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.07_fifo1_rbuf8_level56_buffered_times.csv`；run_id：`20260819_163337`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7619 | 2 | 1.01 | 7575 | 44 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 93 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 733 | 1 | 1 | 733 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 4 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8445 个已退休 batch。

## 18. Eb/N0=3.07 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.07_fifo1_rbuf16_level56_buffered_times.csv`；run_id：`20260819_163656`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7628 | 2 | 1.01 | 7574 | 54 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 93 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 720 | 1 | 1 | 720 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 8 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8441 个已退休 batch。

## 19. Eb/N0=3.07 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.07_fifo1_rbuf32_level56_buffered_times.csv`；run_id：`20260819_164023`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7631 | 2 | 1.01 | 7567 | 64 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 93 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 709 | 1 | 1 | 709 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 16 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8433 个已退休 batch。

## 20. Eb/N0=3.07 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.07_rbuf64_level56_buffered_times.csv`；run_id：`20260819_224406`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7632 | 2 | 1.01 | 7552 | 80 | 0 | 0 | 0.02 | 0 |
| FullEarlyStop | 93 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 692 | 1 | 1 | 692 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 32 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8417 个已退休 batch。

## 21. Eb/N0=3.07 dB; FIFO=ON; R_buf=80; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.07_rbuf80_level56_buffered_times.csv`；run_id：`20260819_230543`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7633 | 2 | 1.01 | 7545 | 88 | 0 | 0 | 0.02 | 0 |
| FullEarlyStop | 93 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 683 | 1 | 1 | 683 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 40 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8409 个已退休 batch。

## 22. Eb/N0=3.07 dB; FIFO=ON; R_buf=96; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.07_rbuf96_level56_buffered_times.csv`；run_id：`20260819_230948`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 7644 | 2 | 1.01 | 7548 | 96 | 0 | 0 | 0.02 | 0 |
| FullEarlyStop | 93 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 664 | 1 | 1 | 664 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 48 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8401 个已退休 batch。

## 23. Eb/N0=3.08 dB; FIFO=ON; R_buf=0; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.08_fifo1_rbuf0_level56_buffered_times.csv`；run_id：`20260819_164348`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8019 | 1 | 1 | 8019 | 0 | 0 | 0 | 0 | 0 |
| FullEarlyStop | 120 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 310 | 1 | 1 | 310 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 24. Eb/N0=3.08 dB; FIFO=ON; R_buf=8; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.08_fifo1_rbuf8_level56_buffered_times.csv`；run_id：`20260819_164723`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8088 | 2 | 1.01 | 8020 | 68 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 120 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 237 | 1 | 1 | 237 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 4 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8445 个已退休 batch。

## 25. Eb/N0=3.08 dB; FIFO=ON; R_buf=16; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.08_fifo1_rbuf16_level56_buffered_times.csv`；run_id：`20260819_165038`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8100 | 2 | 1.01 | 8017 | 83 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 120 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 221 | 1 | 1 | 221 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 8 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8441 个已退休 batch。

## 26. Eb/N0=3.08 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.08_fifo1_rbuf32_level56_buffered_times.csv`；run_id：`20260819_165349`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8101 | 2 | 1.01 | 8010 | 91 | 0 | 0 | 0.02 | 0 |
| FullEarlyStop | 120 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 212 | 1 | 1 | 212 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 16 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8433 个已退休 batch。

## 27. Eb/N0=3.08 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.08_rbuf64_level56_buffered_times.csv`；run_id：`20260819_231259`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8103 | 2 | 1.01 | 7996 | 107 | 0 | 0 | 0.02 | 0 |
| FullEarlyStop | 120 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 194 | 1 | 1 | 194 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 32 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8417 个已退休 batch。

## 28. Eb/N0=3.09 dB; FIFO=ON; R_buf=64; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_rbuf64_level56_buffered_times.csv`；run_id：`20260819_231647`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8254 | 2 | 1.01 | 8162 | 92 | 0 | 0 | 0.02 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 45 | 1 | 1 | 45 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 29. Eb/N0=3.09 dB; FIFO=ON; R_buf=80; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_rbuf80_level56_buffered_times.csv`；run_id：`20260819_232033`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8263 | 2 | 1.01 | 8163 | 100 | 0 | 0 | 0.02 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 36 | 1 | 1 | 36 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 30. Eb/N0=3.09 dB; FIFO=ON; R_buf=128; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_rbuf128_level56_buffered_times.csv`；run_id：`20260819_232417`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8275 | 2 | 1.01 | 8151 | 124 | 0 | 0 | 0.03 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 5 | 1 | 1 | 5 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 19 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8430 个已退休 batch。

## 31. Eb/N0=3.09 dB; FIFO=ON; R_buf=136; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_rbuf136_level56_buffered_times.csv`；run_id：`20260819_233726`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8275 | 2 | 1.02 | 8147 | 128 | 0 | 0 | 0.03 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 1 | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |

> 注意：该运行有 23 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8426 个已退休 batch。

## 32. Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_rbuf138_level56_buffered_times.csv`；run_id：`20260824_172806`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8275 | 2 | 1.02 | 8146 | 129 | 0 | 0 | 0.03 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 注意：该运行有 24 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8425 个已退休 batch。

## 33. Eb/N0=3.09 dB; FIFO=ON; R_buf=140; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_rbuf140_level56_buffered_times.csv`；run_id：`20260819_233054`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8275 | 2 | 1.02 | 8146 | 129 | 0 | 0 | 0.03 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 注意：该运行有 24 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8425 个已退休 batch。

## 34. Eb/N0=3.09 dB; FIFO=ON; R_buf=144; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_rbuf144_level56_buffered_times.csv`；run_id：`20260819_232736`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8275 | 2 | 1.02 | 8146 | 129 | 0 | 0 | 0.03 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 注意：该运行有 24 个 arrival batch 在观测结束时尚未退休（右截尾）；上表的服务次数分布只包含 8425 个已退休 batch。

## 35. Eb/N0=3.1 dB; FIFO=ON; R_buf=32; HISO/SISO=8/8; Group4; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.1_rbuf32_level56_buffered_times.csv`；run_id：`20260819_200826`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8243 | 2 | 1.01 | 8199 | 44 | 0 | 0 | 0.01 | 0 |
| FullEarlyStop | 189 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 17 | 1 | 1 | 17 | 0 | 0 | 0 | 0 | 0 |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。

## 36. Eb/N0=3.09 dB; FIFO=ON; R_buf=138; HISO/SISO=32/32; GlobalPriority; entries=default; drain=OFF

数据源：`/home/zsr71/projects/newcode_level56_buffered_fifo/data/level56_schedule/debug_L6_ebn03.09_fifo1_rbuf138_hiso32_siso32_schedglobal_eqobs1_level56_buffered_times.csv`；run_id：`20260824_194542`。

| 退休原因 | batch数 | max服务 | mean服务 | 服务1次 | 服务2次 | 服务3次 | 服务4次以上 | mean max_pending | P95 max_pending |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Normal | 8299 | 1 | 1 | 8299 | 0 | 0 | 0 | 0 | 0 |
| FullEarlyStop | 150 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| ForcedEvicted | 0 | — | — | 0 | 0 | 0 | 0 | — | — |

> 该运行的 arrival batch 全部已退休；服务次数分布没有尾部截尾。
