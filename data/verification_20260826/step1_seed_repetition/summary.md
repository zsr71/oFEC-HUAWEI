# 第一步第 1a 阶段：HISO/SISO 多 seed 配对复现结果

日期：2026-08-26  
实验目标：在无 FIFO forced eviction 的条件下，检验 `0HISO+8SISO` 与 `8HISO+8SISO` 的 BER 差异是否可在独立随机样本中重复。

## 固定条件

```text
Eb/N0 = 3.065 dB
OFEC_NUM_INFO_BITS = 60,014,592
每 run 的 BER compared bits = 58,608,000
固定 BER cut = 1,406,592 bit
Level 5/6 shared = ON
Group4 scheduler = ON，Group4 entry budget = 8
FIFO = ON，R_buf = 16,000 block row，frame-end drain = ON
SISO = 8；对照变量为 HISO = 0 或 8
```

每个 seed 的两侧运行共用相同的 bit-generator seed 和 channel seed。BER 采用同一全帧固定裁剪口径；没有采用仅统计“已完成区域”的额外筛选。

## 配对结果

| seed | bit/channel seed | 0HISO+8SISO errors | 8HISO+8SISO errors | 0HISO+8SISO BER | 8HISO+8SISO BER | errors ratio | forced eviction（两侧） |
|---|---:|---:|---:|---:|---:|---:|---:|
| seed01 | 20260319 / 3182026 | 11 | 159 | 1.87688e-7 | 2.71294e-6 | 14.45× | 0 / 0 |
| seed02 | 731941 / 810271 | 3 | 77 | 5.11876e-8 | 1.31381e-6 | 25.67× | 0 / 0 |
| seed03 | 1357911 / 2468021 | 19 | 136 | 3.24188e-7 | 2.32050e-6 | 7.16× | 0 / 0 |
| seed04 | 987653 / 123457 | 46 | 375 | 7.84876e-7 | 6.39844e-6 | 8.15× | 0 / 0 |
| seed05 | 440021 / 770003 | 35 | 314 | 5.97188e-7 | 5.35763e-6 | 8.97× | 0 / 0 |
| **合计** | — | **114** | **1061** | **3.89025e-7** | **3.62067e-6** | **9.31×** | **所有 10 runs 为 0** |

合计相比 bit 数为 `293,040,000`（`5 × 58,608,000`）。

## 最终错误位置的配对关系

| seed | 两侧共同错误 | 仅 8HISO+8SISO | 仅 0HISO+8SISO |
|---|---:|---:|---:|
| seed01 | 6 | 153 | 5 |
| seed02 | 1 | 76 | 2 |
| seed03 | 8 | 128 | 11 |
| seed04 | 33 | 342 | 13 |
| seed05 | 22 | 292 | 13 |
| **合计** | **70** | **991** | **44** |

混合组相对纯 SISO 组的净新增错误位置数为 `991 − 44 = 947`。这与按错误计数得到的差值 `1061 − 114 = 947` 一致。

## 本阶段结论与边界

已证实：在本页固定的 3.065 dB、5 对独立 seed 和当前实现配置下，`8HISO+8SISO` 的最终 Post-FEC BER 在每一对中均差于 `0HISO+8SISO`。所有 run 均 `forced_evicted=0`，因此 FIFO 的强制顶出不会解释这个差异。

尚未证实：这并不能证明所有 Eb/N0、所有 HISO 实现或任意 hybrid 类别都有同样结果；也不能证明 `OneMain` 是根因。下一步应按 hybrid 分类禁止 HISO 准入进行受控消融，随后再跟踪首次 writeback/history/后续输入分叉。

## 可复查原始文件

本目录中每个 `*_stdout.log` 保存完整 run 输出；同名 `*_post_fec_errors.txt` 保存该 run 的最终 Post-FEC 错误物理位置。文件名中编码了 seed、HISO 数和 SISO 数。
