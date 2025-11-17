# 接收端：QPSK/BPSK 解调与 LLR 生成

本文说明 `run_pipeline()` 接收端中从 QPSK/BPSK 符号恢复 LLR 的全过程。

## 噪声建模与方差换算

调制符号经 AWGN 信道：  
$$
y_k = x_k + n_k,\quad n_k \sim \mathcal{CN}(0, \sigma^2)
$$
其中 `qam_modulate()` 已将星座归一化到单位平均能量 \(E_s = 1\)。噪声标准差由 `ebn0_to_sigma()` 计算：
$$
\sigma = \sqrt{\frac{N_0}{2}},\quad
N_0 = \frac{1}{E_b / N_0} ,\quad
$$


## BPSK 情况

当 `n_bps = 1` 时，`qam_llr_logsumexp()` 走 BPSK 快捷路径：
$$
L_k = \frac{2}{\sigma^2} \Re\{y_k\}
$$

## QPSK（及偶数位 M-QAM）情况

对于比特数为偶数的调制（QPSK、16QAM、64QAM 等），`build_constellation()` 按 Gray 映射构造星座集合 $\mathcal{S}$。  
对第 $k$ 个接收符号 $y_k$，其第 $b$ 个比特位的 LLR 使用 log-sum-exp 形式：

$$
L_{k,b} =
\log \frac{
  \displaystyle\sum_{s \in \mathcal{S}_b^{(0)}} \exp\!\left(-\frac{\lVert y_k - s \rVert^2}{N_0}\right)
}{
  \displaystyle\sum_{s \in \mathcal{S}_b^{(1)}} \exp\!\left(-\frac{\lVert y_k - s \rVert^2}{N_0}\right)
},
$$

其中 $\mathcal{S}_b^{(0)}$ 表示二进制展开的第 $b$ 位为 0 的所有星座点子集，$\mathcal{S}_b^{(1)}$ 同理为第 $b$ 位为 1 的子集。

在实现中，代码将 `log-sum-exp` 写成“两次遍历 + 最大值偏移”以避免指数下溢：

1. **第一次遍历**：在 $\mathcal{S}_b^{(0)}$ 与 $\mathcal{S}_b^{(1)}$ 中分别找到最大度量
   $$
   m_0 = \max_{s \in \mathcal{S}_b^{(0)}} m_s,\quad
   m_1 = \max_{s \in \mathcal{S}_b^{(1)}} m_s,
   $$
   其中 $m_s = -\lVert y_k - s \rVert^2 / N_0$。
2. **第二次遍历**：在最大值附近累加“归一化后的指数和”
   $$
   \text{sum0} = \sum_{s \in \mathcal{S}_b^{(0)}} \exp(m_s - m_0),\quad
   \text{sum1} = \sum_{s \in \mathcal{S}_b^{(1)}} \exp(m_s - m_1).
   $$
3. **合成 LLR**：
   $$
   L_{k,b} = \bigl(m_0 + \log \text{sum0}\bigr) - \bigl(m_1 + \log \text{sum1}\bigr).
   $$

## 实现细节

- `qam_llr_from_ebn0()` 先由给定的 $E_b/N_0$ 计算噪声标准差 `sigma`，并复用同一套公式驱动 BPSK、QPSK 以及更高阶偶数位 QAM。
- 随后 `apply_known_zero_prefix()` 会基于 LLR 的符号和幅度施加先验约束，为 oFEC 解码提供一致的输入软信息。

通过上述步骤，发射端使用的 Gray 映射与接收端 LLR 推导在数学上严格匹配，从而保证 BPSK/QPSK 及更高阶偶数位 QAM 仿真在不同 $E_b/N_0$ 下具有良好的可重复性。
