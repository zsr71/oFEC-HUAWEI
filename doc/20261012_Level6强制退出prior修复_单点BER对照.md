本次在当前双窗口工程上实施 Level 6 最终强制退出 prior 修复，并运行完整编码、BPSK、AWGN、量化、译码、信息提取的修复前后成对实验。工作点的 Post-FEC BER 从 1.3574938575×10⁻⁴ 降至 4.5300982801×10⁻⁵，净减少 5,301 个错误，下降 66.628959%。

## 实验配置

| 项目 | 设置 |
|---|---|
| 工程 | newcode_level56_dual_window |
| 运行模式 | Level 5/6 双窗口独立浮动 |
| Eb/N0 | 3.06 dB |
| R_buf | 32 |
| HISO/SISO | 8/8 |
| LLR | 6 bit，clip ratio 0.5 |
| 核心 | chase_baseline，CHASE_L=6，CHASE_NTEST=64 |
| 调度 | Group4LoadSortedMultiround，8 entries，Level5First，AllGroups |
| 信息/噪声种子 | 20260319 / 3182026 |
| 信息长度 | 60,014,592 bit，256 个完整窗口 |
| 编码长度 | 69,476,352 bit |
| BER 统计范围 | 固定裁掉前 4 个和后 2 个完整窗口；58,608,000 bit |
| 帧尾 | dual drain=1；FIFO drain=0 |
| 编译 | GCC 11.4，CMake Release |

before 是本次实验开始前的完整源码快照，包含工程原有未提交修改；after 只叠加本次修复。源码快照在 `/tmp/level56-ber-before-20261012`。二者重新编译运行，修复前 7,956 错误也精确复现了原长帧记录。

## 结果

| 指标 | 修复前 | 修复后 |
|---|---:|---:|
| 比较信息比特 | 58,608,000 | 58,608,000 |
| Pre-FEC 错误 | 1,295,726 | 1,295,726 |
| 量化硬判错误 | 1,298,985 | 1,298,985 |
| Post-FEC 错误 | 7,956 | 2,655 |
| Post-FEC BER | 0.000135749385749 | 0.000045300982801 |

BER 绝对下降 9.04484029484×10⁻⁵，错误数约减少到原来的 1/3（比值 2.99661）。

逐位交集比较：原来错误、修复后正确 6,043 bit；原来正确、修复后错误 742 bit；两次都错误 1,913 bit。因此 before=6,043+1,913=7,956，after=742+1,913=2,655，净收益=6,043−742=5,301。

补存 prior 并不保证逐位单调改善：prior 本身可能错误或幅度足以翻转原本正确的信道硬判。本次仍有 742 个新增错误，净改善来自消除的 6,043 个错误更多。

## 修复和归因验证

1. Pending 仍可服务时不提交最终 history。
2. AlreadyDecoded 且 writeback_complete 时保留已经产生的最终 history。
3. Level 6 最终 ForcedEvicted 且未完成写回时，按 writeback 的历史半码字 128-bit 地址映射，从当前 work_llr 复制 prior 到最终 history；不乘 alpha，不加入新的 extrinsic。最终输出保持原公式 channel+history。

修复挂在 dual boundary 的 Level 6 force 路径；帧尾 drain 共用该路径。同工程的原绑定 FIFO retirement 路径也修复，但本次未单独运行其 BER 对照。其他五个旧工程没有在本次实验中传播修复，不能把本结果当作它们的实测结果。

新增回归涵盖 pending 不提交、完成结果不覆盖、Level 5 不提交 Level 6 history、只写 128 个历史映射位置、float/qfloat<6> 转换及双窗口/FIFO 强制退出入口。全部现有共享调度回归也通过。

逐行比较 CSV（仅排除 run_id 和 label）：1,983,552 条 code 记录、34,165 条共享调度轮次、70,333 条双窗口微步、16,935 条服务时刻记录，差异全部为 0。因此本次没有增加 HISO/SISO 资源或改变译码调度。默认关闭逐 code 输入/输出 hash 观测；本次没有以未导出的 hash 作证据。

Level 6 强制退出 625 批，退出时有 2,640 个未完成写回码字。历史地址映射在 BER 区间内覆盖 279,520 个信息位置。全部 6,043 个被纠正位置和全部 742 个新增错误位置均位于这个映射集合中。

实际位置例子：信息比特索引 3019208，对应 B905 的 Level 6 码字全局行 29656、BCH bit k=8，其 history 地址是 (29312, 8)。这个位置修复前错误，补存当前 prior 后正确。新增错误例子：信息比特 3149674，对应 B938 的码字行 30689、k=54，history 地址 (30487, 49)。这些地址根据实际退出与调度记录重建，未采集该位置的数值 LLR，不能据此声称具体幅度。

## 复现和原始文件

在工程根目录运行：

```bash
./data/prior_fix_ber_20261012/run_pair.sh before dual
./data/prior_fix_ber_20261012/run_pair.sh after dual
python3 data/prior_fix_ber_20261012/analyze_pair.py dual
python3 data/prior_fix_ber_20261012/attribute_forced_positions.py
```

原始日志：`data/prior_fix_ber_20261012/dual_before/stdout.log` 与 `dual_after/stdout.log`。对应目录内有 timing.log、errors.csv 和调度 CSV。汇总见 dual_results.csv/json；位置明细见 dual_changed_position_examples.csv。experiment_only.patch 仅含本实验改动；full_source_manifest.json 保存全部源码和二进制 SHA-256，regression.log 保存回归结果。

修复前/后运行时间分别为 5:55.83 / 5:45.27，峰值 RSS 均约 5.48 GB。两者运行部分重叠，耗时不适合作性能改善结论。

这是一个 Eb/N0、一个种子、一个配置上的成对样本结果，足以量化该样本的修复影响；不是多种子 BER 估计，也不能把 66.63% 外推到所有工作点或视作 dB 编码增益。错误有码字/帧聚集，因此不套用独立 bit 假设给出过窄置信区间。
