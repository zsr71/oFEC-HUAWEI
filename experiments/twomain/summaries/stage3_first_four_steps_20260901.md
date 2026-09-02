# TwoMain 阶段三前四步实施与验证记录

日期：2026年9月1日

本记录只覆盖阶段三的前四步：实现方案六统一参数化固定输出、保持默认 Legacy 行为、复核 `B134 / Level 5 / code 4`、完成短帧冒烟。尚未开始正式的 $M_2$ 长帧扫描，也未开始方案四、方案五。

## 一、实现结果

方案六-A、六-B和六-C使用同一条公式和同一条代码路径：

\[
M_i=
\begin{cases}
\rho_{\mathrm{corr}}M_2, & i\in\mathcal E,\\
\rho_{\mathrm{keep}}M_2, & i\notin\mathcal E,
\end{cases}
\]

\[
L_{\mathrm{post},i}=s_iM_i,
\]

\[
L_{\mathrm{out},i}=L_{\mathrm{post},i}-L_{\mathrm{in},i}.
\]

其中：

1. $\mathcal E$ 是 BCH 译码器相对输入硬判实际翻转的两个主体位置，不依赖真实发送码字；
2. $s_i$ 由 BCH 纠正后的合法硬码字决定；
3. 六-A取 $\rho_{\mathrm{corr}}=\rho_{\mathrm{keep}}=1$；
4. 当时定义的六-B方向一固定 $M_2$，降低 $\rho_{\mathrm{keep}}$；2026-09-02另增固定$\rho_{\mathrm{keep}}=1$、降低$\rho_{\mathrm{corr}}$的方向二；
5. 六-C继续联合调整两个系数；
6. 当前公式中系数为零的精确定义是零后验，即 $L_{\mathrm{out}}=-L_{\mathrm{in}}$。本轮定点检查和冒烟没有使用零系数。

新增公共参数为：

```text
TWOMAIN_HISO_OUTPUT_MODE
TWOMAIN_HISO_M2
TWOMAIN_HISO_RHO_CORR
TWOMAIN_HISO_RHO_KEEP
```

单点程序提供对应运行时入口：

```text
OFEC_TWOMAIN_HISO_OUTPUT_MODE=legacy|parameterized
OFEC_TWOMAIN_HISO_M2=<非负有限数>
OFEC_TWOMAIN_HISO_RHO_CORR=<0到1>
OFEC_TWOMAIN_HISO_RHO_KEEP=<0到1>
```

参数化输出只在以下条件同时成立时生效：

```text
类别为 TwoMain
最终动作是 HISO
输出模式为 parameterized
```

其他类别继续调用原来的 HISO 输出函数。EarlyStop、分类、HISO准入、MUX、资源调度、BCH纠正、合法性检查、Level 5/6写回、Level 6 history、FIFO和pending逻辑均未修改。

## 二、Legacy 默认行为验证

默认模式显式设为 `Legacy`，并显式进入原输出函数，不依靠 $M_2=99$、两个系数为1来近似旧行为。

已通过两层回归：

1. `ofec_level56_shared_regression_check` 通过；
2. `check_twomain_default_equivalence.py --allow-dirty` 通过。

冻结签名仍为：

| 观测量 | 结果 |
|---|---:|
| Post-FEC | `0 / 412032` |
| 完成 batch | `513` |
| forced eviction | `0` |
| 最大 FIFO 深度 | `18` |
| Level 5动作 | `[13765, 598, 2053, 576]` |
| Level 6动作 | `[16361, 0, 55, 576]` |
| entry/HISO/SISO/idle总量 | `[2108, 598, 2108, 3650, 2140]` |

Post-FEC错误位置文件仍为空，哈希仍为：

```text
e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855
```

这说明只增加参数能力、但不启用方案六时，现有解码结果没有发生变化。

## 三、B134定点检查

### 3.1 复现实验条件

```text
Eb/N0 = 3.065 dB
NUM_INFO_BITS = 60014592
bitgen seed = 20260319
channel seed = 3182026
FIFO = ON
R_buf = 16000
frame-end drain = ON
HISO类别掩码 = 15
HISO/SISO = 8/8
Group4 entry = 8
目标 = B134 / Level 5 / code 4
```

不能用短帧中的同编号 `B134` 代替该样本。短帧会改变整帧动态量化裁剪值，因而相同编号的输入状态不同；本轮已经按原6000万bit条件重新复现。

### 3.2 Legacy复现结果

Legacy定点结果与历史记录一致：

```text
buffered batch = B134
source level = 5
batch code = 4
hybrid class = TwoMain
final action = HISO
输入硬判相对真实码字错误数 = 2
BCH实际翻转位置 = 153,199
BCH纠正后相对真实码字错误数 = 0
```

新增参数说明行之外，Legacy trace与历史trace一致。完整长帧结果也再次复现为：

```text
Post-FEC = 159 / 58608000
forced eviction = 0
最大 FIFO 深度 = 1888
```

### 3.3 参数化模式局部结果

局部检查使用：

```text
M2 = 32
rho_corr = 1
rho_keep = 0.5
```

与Legacy相比，下列内容保持一致：

1. `lch`输入一致；
2. `lin`输入一致；
3. 类别仍为TwoMain；
4. 动作仍为HISO；
5. BCH实际翻转位置仍为 `153,199`；
6. BCH纠正后硬码字完全一致，且相对真实码字为0错；
7. 写回目标坐标不变。

第一处变化发生在 `lout` 生成处。两个纠正位置经过原有后处理和量化后的观测为：

| 位置 | 写回坐标 | Legacy `lout` | 参数化 `lout` |
|---:|---|---:|---:|
| 153 | `(4963,26)` | `-1.00253522` | `-0.711476564` |
| 199 | `(4963,68)` | `+1.00253522` | `+0.614457071` |

因此，本轮已经确认方案六没有改变TwoMain之前的分类、调度或BCH纠正，其第一处分叉确实是TwoMain外信息输出。

需要强调：同一组参数运行完整长帧时出现严重劣化，Post-FEC达到 `3190321 / 58608000`，并产生6999个forced-evicted batch。该运行只证明参数化输出会沿后续迭代传播并真实生效；它同时表明 `M2=32、rho_keep=0.5` 不是可接受候选，不能据此评价方案六整体。

## 四、短帧冒烟结果

三组均使用固定短帧、相同信道和随机种子：

| 模式 | 参数 | 正常退出 | Post-FEC | forced eviction | 最大 FIFO 深度 |
|---|---|---|---:|---:|---:|
| Legacy | 原始固定输出 | 是 | `0 / 412032` | 0 | 18 |
| 六-A | $M_2=32$，两个系数均为1 | 是 | `0 / 412032` | 0 | 30 |
| 六-B机制冒烟 | $M_2=32$，$\rho_{\mathrm{corr}}=1$，$\rho_{\mathrm{keep}}=0.5$ | 是 | `44443 / 412032` | 47 | 257 |

短帧结果的用途是检查程序能否正常运行、参数是否生效、默认是否不变，不用于筛选最终BER参数。六-B这组激进参数已经在短帧中表现出明显劣化，后续不进入正式候选集合。

对应可复现实验目录为：

```text
runs/twomain/regression_default_quick/stage3_smoke_legacy_20260901
runs/twomain/scheme6_uniform_32_smoke/stage3_smoke_6a_20260901
runs/twomain/scheme6_differential_32_smoke/stage3_smoke_6b_20260901
```

## 五、本轮新增验证能力

快速回归现在会检查：

1. Legacy TwoMain输出逐位置保持旧公式；
2. 两错BCH能够恢复合法原码字；
3. 六-A对全部位置使用统一幅度；
4. 六-B和六-C能区分两个实际纠正位置与其他位置；
5. 纠正位置集合恰好包含两个主体位置；
6. 方案六参数不影响非TwoMain类别；
7. $M_2$和两个系数的运行时范围受到校验。

实验运行器已经能够从受版本控制的配置文件传递方案六参数。本轮增加：

```text
experiments/twomain/specs/scheme6_uniform_32_smoke.toml
experiments/twomain/specs/scheme6_differential_32_smoke.toml
```

## 六、剩余工作

阶段三下一步从第五步开始：只做六-A的单变量 $M_2$ 长帧扫描。尚未执行的内容是：

1. 选择一组不过度激进的 $M_2$ 序列；
2. 在固定长帧和固定随机种子下运行六-A；
3. 同时比较BER、错误位置、TwoMain动作数、早停率、FIFO深度和forced eviction；
4. 只保留少量有希望的 $M_2$；
5. 之后再以相同 $M_2$ 比较六-B；
6. 只有六-B显示独立收益时才进入六-C；
7. 最后才进行多随机种子和多Eb/N0验证。

本轮没有提交或推送Git，也没有删除任何既有运行数据。
