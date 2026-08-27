# Eb/N0=3.1 dB、R_buf=32：五六级 Buffered FIFO 运行分析

## 演示文稿

`level56_buffered_fifo_ebn03p1_rbuf32_analysis.pptx`

共 13 页，沿用项目已有的 16:9、白底、红/蓝强调色演示风格。

## 数据来源

本次运行命令：

```bash
OFEC_EBN0_DB=3.1 LEVEL56_BUFFER_ROWS=32 ./build/ofec_single
```

对应输出：

- `data/run_20260819-200826_single.log`
- `data/level56_schedule/debug_L6_ebn03.1_rbuf32_level56_buffered_times.csv`
- `data/level56_schedule/ofec_single_ebn03.1_rbuf32_level56_schedule_rounds.csv`
- `data/level56_schedule/ofec_single_ebn03.1_rbuf32_level56_schedule_codes.csv`

## 主要实测数字

| 指标 | 数值 |
|---|---:|
| Eb/N0 | 3.1 dB |
| R_buf | 32 block row |
| Post-FEC BER | 2.09785e-7 |
| Post-FEC errors | 6 / 28,600,704 |
| 服务时刻 | 8,449 |
| C_t=0 | 61 |
| C_t=1 | 8,346 |
| C_t>=2 | 42 |
| 全 EarlyStop batch | 189 |
| ForcedEvicted batch | 17 |
| 最大 FIFO depth | 17（服务前，单位为 batch） |

## R_buf=32 与 R_buf=64 对照

在相同 `Eb/N0=3.1 dB`、种子、共享资源和调度下，`R_buf=64` 的单点结果为：

| 指标 | R_buf=32 | R_buf=64 |
|---|---:|---:|
| ForcedEvicted batch | 17 | 0 |
| 最小 S_t | 0 | 6 |
| 最大 FIFO depth（服务前） | 17 | 30 |
| 最长正常等待 | 16 t | 29 t |
| Post-FEC errors | 6 | 5 |
| Post-FEC BER | 2.09785e-7 | 1.74821e-7 |

`R_buf=32` 下被强制退休的 17 个 batch，在 `R_buf=64` 下都变为 `Normal` 完成。典型例子：

- B2985：R_buf=32 在 `t=3001` ForcedEvicted；R_buf=64 在 `t=3002` Normal 完成。
- B3535：R_buf=32 在 `t=3551` ForcedEvicted；R_buf=64 在 `t=3564` Normal 完成，等待 29 t，仍保留 6 行 buffer 余量。

## 讲解重点

建议重点讲第 4、5、7、8、10、11 页：

1. 第 4 页说明每个 t 固定输出两行，C_t 只决定窗口如何运动。
2. 第 5 页区分 FIFO depth、队首 batch 和物理输出行数。
3. 第 7 页说明 C_t=0 时队首不变，只有完成后才切换。
4. 第 8 页说明 S_t=0 时 ForcedEvicted 与固定两行输出节拍的关系。
5. 第 10 页说明 R_buf 翻倍后，17 个 ForcedEvicted 全部消失。
6. 第 11 页用 B2985 / B3535 展示同一 batch 在 32 与 64 行 buffer 下的不同命运。

## 术语边界

- `FIFO depth` 是 FIFO 中 64-code batch 的数量，不是 block row 数。
- `forced_evicted_global_rows` 是被强制退休的 batch 中仍未完成 code 的诊断行列表，不是物理 SRAM 顶出行数。
- 物理 SRAM 每个时刻仍固定输出两行；ForcedEvicted 是同一时刻发生的逻辑 batch 退休事件。
- 本次状态序列可以验证 FIFO、C_t、S_t 和 ForcedEvicted 的逻辑；当前软件仿真仍需单独核对连续 SRAM 地址搬移是否完整接入读写路径。
