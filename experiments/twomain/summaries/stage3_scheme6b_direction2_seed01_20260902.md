# 阶段三：方案六-B方向二固定长帧seed01筛选

## 一、实验目的

既有方案六-B只测试了方向一：保持两个BCH纠正位置的后验幅度，降低其余254个非纠正位置的后验幅度。方向一在$M_2=48$时发生闭环失稳，不能据此否定相反的位置差异化。

本轮测试方向二：保持254个非纠正位置的幅度，只降低两个BCH实际翻转位置的幅度。目的是回答：

> 在保持系统整体稳定的前提下，减弱两个纠正位置的后验强度，能否进一步改善六-A统一幅度的BER？

## 二、公式和参数

沿用方案六统一公式：

$$
M_i=
\begin{cases}
\rho_{\mathrm{corr}}M_2,&i\in\mathcal E,\\
\rho_{\mathrm{keep}}M_2,&i\notin\mathcal E,
\end{cases}
$$

$$
L_{\mathrm{out},i}=s_iM_i-L_{\mathrm{in},i}.
$$

固定：

$$
M_2=48,
\qquad
\rho_{\mathrm{keep}}=1,
$$

扫描：

$$
\rho_{\mathrm{corr}}\in\{0.75,0.5,0.25\}.
$$

实验与六-A `1/1`采用相同的固定长帧seed01：

```text
Eb/N0 = 3.065 dB
bitgen/channel seed = 20260319 / 3182026
信息bit数 = 60,014,592
BER比较bit数 = 58,608,000
Level 5/6 shared = ON
HISO/SISO = 8/8
Group4 entry = 8
FIFO = ON
R_buf = 16,000 block row
frame-end drain = ON
```

除实验身份和`rho_corr`外，三份方向二配置与六-A `M2=48`配置完全一致。

## 三、BER和FIFO结果

| 配置 | $M_{\mathrm{corr}}$ | $M_{\mathrm{keep}}$ | Post-FEC错误 | BER | forced eviction | 最大FIFO深度 | 服务时刻数 | 最小$S_t$ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 六-A：`1/1` | 48 | 48 | 64 | $1.09200\times10^{-6}$ | 0 | 1869 | 19082 | 11630 |
| 方向二：`0.75/1` | 36 | 48 | 85 | $1.45031\times10^{-6}$ | 0 | 1876 | 19095 | 11604 |
| 方向二：`0.5/1` | 24 | 48 | 74 | $1.26263\times10^{-6}$ | 0 | 1890 | 19130 | 11534 |
| 方向二：`0.25/1` | 12 | 48 | 75 | $1.27969\times10^{-6}$ | 0 | 1911 | 19189 | 11416 |

三个方向二点都正常完成16897个batch，且均无forced eviction。这与方向一`M_2=48, 1/0.5`的3765次forced eviction形成明确对照，说明旧方向一的失稳主要与254个keep位置的有效幅度降到24有关，而不是“只要纠正位和非纠正位不同就会失稳”。

但是，三个方向二点均未优于六-A统一幅度的64错：

```text
rho_corr=0.75：85错，比六-A多21错
rho_corr=0.50：74错，比六-A多10错
rho_corr=0.25：75错，比六-A多11错
```

方向二内部不是严格单调关系。`0.5/1`是本轮最好点，但仍未达到进入五组配对随机种子的门槛。

## 四、错误位置差集

相对六-A `1/1`的64个最终错误位置：

| 方向二配置 | 共同错误 | 仅方向二错误 | 仅六-A错误 |
|---|---:|---:|---:|
| `0.75/1` | 37 | 48 | 27 |
| `0.5/1` | 33 | 41 | 31 |
| `0.25/1` | 26 | 49 | 38 |

方向二同样不是简单地保留或删除六-A错误集合。改变两个corr位置的外信息后，后续EarlyStop、HISO/SISO动作和history输入发生闭环传播，最终错误集合重新分布。

## 五、TwoMain动作和调度

| 配置 | TwoMain HISO | TwoMain SISO | TwoMain pending | Level 5 Unscheduled | Level 6 Unscheduled |
|---|---:|---:|---:|---:|---:|
| 六-A：`1/1` | 2789 | 13516 | 8370 | 72128 | 72128 |
| 方向二：`0.75/1` | 2800 | 13525 | 8384 | 72544 | 72544 |
| 方向二：`0.5/1` | 2828 | 13603 | 8645 | 73664 | 73664 |
| 方向二：`0.25/1` | 2919 | 13733 | 9073 | 75520 | 75520 |

随着`rho_corr`降低，服务时刻数、最大FIFO深度和pending总体缓慢增加，但没有进入失稳区。动作变化是外信息改变后的闭环结果，并非配置直接改变了TwoMain准入或调度规则。

## 六、结论和后续决定

本轮能够确定：

1. 保持$\rho_{\mathrm{keep}}=1$能够避免旧方向一的大规模FIFO失稳；
2. 因此，旧方向一的失败不能解释成“所有位置差异化都会失稳”；
3. 但在$M_2=48$和seed01下，降低$\rho_{\mathrm{corr}}$没有获得超过六-A `1/1`的BER收益；
4. 三个方向二点均不进入五组配对随机种子；
5. 不依据本轮结果重启六-C；
6. 方案六当前最佳配置仍为六-A：$M_2=48,\rho_{\mathrm{corr}}=\rho_{\mathrm{keep}}=1$；
7. 方案六-B两个方向均已完成预定的第一轮验证，固定位置差异化暂未显示相对六-A的额外收益。

本轮结果不证明任意连续取值的`rho_corr`都不可能优于六-A；它证明的是预先规定的三个代表点均未达到推进门槛。基于当前研究顺序，不继续对单seed做更密集的参数拟合，以避免围绕单一随机样本过拟合。

## 七、配置和运行位置

配置：

```text
experiments/twomain/specs/scheme6b_dir2_m2_48_corr75_keep100_long.toml
experiments/twomain/specs/scheme6b_dir2_m2_48_corr50_keep100_long.toml
experiments/twomain/specs/scheme6b_dir2_m2_48_corr25_keep100_long.toml
```

运行：

```text
runs/twomain/scheme6b_dir2_m2_48_corr75_keep100_long/stage3_dir2_seed01_20260902/
runs/twomain/scheme6b_dir2_m2_48_corr50_keep100_long/stage3_dir2_seed01_20260902/
runs/twomain/scheme6b_dir2_m2_48_corr25_keep100_long/stage3_dir2_seed01_20260902/
```

本轮运行使用`--allow-dirty`，运行清单记录了完整工作树、参数和可执行文件哈希；可用于阶段内候选筛选，最终正式引用仍应在代码提交后的干净工作区复核。
