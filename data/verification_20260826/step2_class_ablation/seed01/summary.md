# 第二步第 2a 阶段：HISO 类别准入阶梯（seed01）

## 条件

```text
Eb/N0 = 3.065 dB
bit/channel seed = 20260319 / 3182026
OFEC_NUM_INFO_BITS = 60,014,592
BER compared bits = 58,608,000；固定 cut = 1,406,592
FIFO = ON，R_buf = 16,000，frame-end drain = ON
Group4 entry budget = 8
```

所有五组均为 `forced_evicted=0`。

`LEVEL56_HISO_CLASS_MASK` 的 bit 顺序为：`ParityOnly / OneMain / OneMainPlusParity / TwoMain`。被禁用的类别不会完成，而是继续成为 SISO/pending 候选。

## BER 阶梯

| 组 | mask | HISO capacity / SISO capacity | 允许走 HISO 的类别 | errors | BER |
|---|---:|---:|---|---:|---:|
| A | 0x0 | 0 / 8 | 无 | 11 | 1.87688e-7 |
| B | 0x8 | 8 / 8 | TwoMain | 355 | 6.05719e-6 |
| C | 0xC | 8 / 8 | TwoMain、OneMainPlusParity | 271 | 4.62394e-6 |
| D | 0xE | 8 / 8 | TwoMain、OneMainPlusParity、OneMain | 191 | 3.25894e-6 |
| E | 0xF | 8 / 8 | 全类别 | 159 | 2.71294e-6 |

E 组精确复现第一步 seed01 的全类别结果（159 errors）。

## 各类别实际动作数（HISO / SISO / Unscheduled）

| 类别 | B: 0x8 | C: 0xC | D: 0xE | E: 0xF |
|---|---:|---:|---:|---:|
| ParityOnly | 0 / 788 / 416 | 0 / 793 / 427 | 0 / 784 / 205 | 318 / 464 / 202 |
| OneMain | 0 / 104415 / 54999 | 0 / 104287 / 54403 | 29233 / 74413 / 26049 | 28999 / 74681 / 26150 |
| OneMainPlusParity | 0 / 226 / 179 | 118 / 107 / 159 | 40 / 182 / 92 | 40 / 186 / 98 |
| TwoMain | 10126 / 6795 / 14614 | 10097 / 6815 / 14467 | 2924 / 13489 / 8691 | 2851 / 13576 / 8704 |
| HardFail（始终 SISO-only） | 0 / 6005 / 8669 | 0 / 5948 / 8409 | 0 / 5571 / 5260 | 0 / 5579 / 5270 |

## 解读边界

不能把 B 相对 A 的 BER 差异简单归因给 TwoMain 本身。关闭 OneMain 的 HISO 准入后，OneMain 会竞争 SISO 机会，进而使 TwoMain 的实际 HISO 动作数从 E 的 2,851 增至 B 的 10,126；其它类别的 SISO/pending 服务也同时改变。

当前只可说：在这一个 seed 的无 forced-eviction 条件下，允许 TwoMain HISO 的 B 组已显著劣于纯 SISO；而加入其它类别后 BER 反而逐步下降。尚不能据此为任何单一类别下因果结论。

下一轮应从 `0xF` 出发一次只关闭一个类别：

```text
0xE：仅关闭 ParityOnly
0xD：仅关闭 OneMain
0xB：仅关闭 OneMainPlusParity
0x7：仅关闭 TwoMain
```

并为每组同时保存实际 HISO/SISO/Unscheduled 动作数与 post-FEC error position 集合。
