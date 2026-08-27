# 第二步第 2b 阶段：全类别基线逐一关闭一类 HISO（seed01）

## 固定条件

```text
Eb/N0 = 3.065 dB
bit/channel seed = 20260319 / 3182026
OFEC_NUM_INFO_BITS = 60,014,592
BER compared bits = 58,608,000；固定 cut = 1,406,592
FIFO = ON，R_buf = 16,000，frame-end drain = ON
Group4 = 8 entries，物理资源固定为 8HISO + 8SISO
```

所有 run 的 `forced_evicted=0`。

类别掩码 bit 顺序：`ParityOnly / OneMain / OneMainPlusParity / TwoMain`。

## BER 结果

| 配置 | mask | 关闭类别 | errors | BER | 与全类别共同 / 仅全类别 / 仅本配置错误位置 |
|---|---:|---|---:|---:|---:|
| 全类别基线 | 0xF | 无 | 159 | 2.71294e-6 | — |
| no_parityonly | 0xE | ParityOnly | 191 | 3.25894e-6 | 128 / 31 / 63 |
| no_onemain | 0xD | OneMain | 308 | 5.25526e-6 | 58 / 101 / 250 |
| no_onemainplusparity | 0xB | OneMainPlusParity | 176 | 3.00300e-6 | 113 / 46 / 63 |
| **no_twomain** | **0x7** | **TwoMain** | **9** | **1.53563e-7** | **4 / 155 / 5** |

## 关键动作数（HISO / SISO / Unscheduled）

| 类别 | 全类别 0xF | 关闭 TwoMain 0x7 |
|---|---:|---:|
| ParityOnly | 318 / 464 / 202 | 323 / 460 / 215 |
| OneMain | 28999 / 74681 / 26150 | 29963 / 73176 / 26344 |
| OneMainPlusParity | 40 / 186 / 98 | 49 / 174 / 99 |
| TwoMain | 2851 / 13576 / 8704 | 0 / 16052 / 8722 |
| HardFail | 0 / 5579 / 5270 | 0 / 5190 / 4891 |

## 结论边界

关闭 TwoMain HISO 后，159 个全类别基线错误中有 155 个消失，仅产生 5 个新错误位置，并且无 forced eviction。由此，TwoMain HISO 是当前最强的高优先级嫌疑类别。

但不能把它写成逐 code 的最终因果证明。原因是关闭 TwoMain 后，其他类别的 HISO/SISO/Unscheduled 数也发生了小幅变化；例如 OneMain HISO 从 28,999 变为 29,963。下一步应对照 `mask=0xF` 和 `mask=0x7` 的同 seed 轨迹，以实际的 `(batch_id, level, code)` action、HISO 输出、tile writeback 和下游 lin 的首次分叉建立因果链。
