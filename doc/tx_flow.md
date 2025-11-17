# 发射端逻辑说明

本文概述仿真中发射端的处理链路，并重点解释 oFEC 矩阵的构成方式。

## 1. 信息位注入与 oFEC 编码

1. `generate_bits(params)` 负责依据 `Params::NUM_INFO_BITS` 和随机种子生成原始信息位流，这些位会被组织成多个 32×111 的“信息矩形”（3552 bit），以满足 oFEC 编码对输入尺寸的对齐要求。
2. `ofec_encode(info_bits, params)` 将信息位编码成一张二维比特矩阵 `Matrix<uint8_t>`。矩阵的列数恒为 `128`，行数根据需要动态增长（初始高度为 `一个tile的高度=22行`，随后一行一行的填入BCH编码后的比特）。

### oFEC 矩阵的组成细节

- **整体坐标系**：矩阵行号等价于 `(R, r)`，即子块行 `R` 与子块内比特行 `r`（`r ∈ [0, 15]`）；列号等价于 `(C, c)`（`C ∈ [0, 7]`，`c ∈ [0, 15]`）。因此一个行向量可视作 8 个 16 bit 子块的拼接。
- **零前缀与窗口对齐**：编码前先构造 `u = [zero_prefix, bits]`。其中 `zero_prefix = PAD_RECTS × 3552`，`PAD_RECTS = 22 / 2 = 11`，保证最初的若干“历史行”是全零，便于 `read_hist_bit()` 读取到已知保护区。
- **左半 128 bit（历史位）**：每生成一行，`read_hist_bit()` 会按照 `{(R^1) − 2G − 2·NB + 2·⌊k/B⌋ , ⌊k/B⌋ , (k mod B) ^ r , r}` 的寻址(参考标准word的第28页)在矩阵中回读 128 个“历史位”。
- **右半 111 bit（系统位）**：`u_index()` 负责将 $k\in[0,110]$ 映射到全局信息比特流 $u$。其映射关系为
$$
\begin{aligned}
W_{R,r}(128 + k)
  &= u\Big(
      \big\lfloor\tfrac{R}{2}\big\rfloor \cdot 32 \cdot 111
      + \big( (R \bmod 2)\cdot 16 + \underline{r} \big)
        \cdot \Bigl(16 - \Big\lfloor\tfrac{k}{96}\Big\rfloor\Bigr) \\
  &\qquad\quad
      + \Big\lfloor\tfrac{k}{16}\Big\rfloor \cdot 512
      + (k \bmod 16)
     \Big),
\end{aligned}
$$
其中 $\lfloor\cdot\rfloor$ 表示向下取整，$\bmod$ 表示取模运算。按照上述规则，每一行准确提取 111 个全新的系统信息位。

- **BCH(255,239) 与校验写回**：
  - 将左 128 bit 与右 111 bit 拼接得到 239 bit 信息字 `msg239`；
  - 通过 `bch_255_239_parity()` 生成 16 个局部校验位，并计算 overall parity；
  - `write_right_to_mat()` 依据和前面相同的寻址，把 111 个系统位（索引 128..238）、16 个 parity（239..254）以及 overall（255）写入当前 `(R, r)` 行中对应列，必要时动态扩行。
- 如此迭代，矩阵逐步写满，为后续调制和 LLR 提供行优先布局。

## 2. 发射侧参考 LLR

编码后立即调用 `hard_bits_to_llr_matrix()` 将 0/1 映射为 ±A（默认 50），再经 `rx_info_from_bit_llr()` 提取 `tx_info_bits_ref`，作为比较 BER 的理想参考。该步骤不参与上行传输，但能保证发射端和接收端使用同一信息提取逻辑，方便后面进行BER的比较。

## 3. 展平、交织与调制

1. `flatten_row_major()` 以行优先顺序把矩阵转成线性比特流 `coded_bits`（`src/pipeline_runner.cpp:134-138`），遵循 oFEC 子块布局。
2. 根据矩阵外形和 `config.interleaver_name` 创建 `Interleaver`，将比特按 16×16 子块粒度做块交织（`src/pipeline_runner.cpp:140-147`），以匹配硬件/链路对时序的要求。
3. 选择调制阶数 `n_bps`（默认 QPSK，可选BPSK），调用 `qam_modulate()` 将交织后比特映射到复符号序列（`src/pipeline_runner.cpp:149-168`）。

## 4. 信道建模

- 使用 `add_awgn()` 按所选 Eb/N0、码率（由 `TAKEBITS / N = 111/128` 推得）和随机种子注入 AWGN，得到接收符号 `rx_syms`（`src/pipeline_runner.cpp:170-185`）。

