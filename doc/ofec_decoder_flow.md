# oFEC 解码整体流程说明

本节总结 `src/ofec_decoder.cpp` 与 `src/decoder_core.cpp` 中的关键步骤，展示窗口化 oFEC 解码的输入、调度及输出组合方式。

## 1. 输入矩阵与先验处理

- 接收端先将串行 LLR 反交织得到 `channel_llr = Matrix<LLR>(rows, 128)`。
- `apply_known_zero_prefix()` 会把 warm-up 保护区及 guard 子行改写成极大正LLR，因为这些行在初始发射的时候就是全0的。
- 若使用量化 LLR，`llr_to_float()` / `llr_from_float()` 始终负责在 `qfloat` 与浮点之间映射。

## 2. 滑动窗口调度

窗口高度定义为
$$
H_\text{win}= T_\text{win} \cdot H_\text{tile}- (T_\text{win} - 1) \cdot \text{overlap},
$$
其中 $T_\text{win}$ 为每个窗口内包含的 tile 数量，$H_\text{tile}$ 为单个 tile 的高度，`overlap` 为相邻 tile 之间的重叠行数。

`process_window_impl()` 从 `win_start` 开始，在输入矩阵中复制出一段高度为 $H_\text{win}$ 的工作区，对应的结束行号为
$$
win_\text{end} = win_\text{start} + H_\text{win} - 1.
$$
在该工作区范围内，函数按窗口内的 tile 顺序逐个调用 `process_tile_impl()` 完成局部译码与外信息更新。


窗口每次滑动 `pop_push_rows = WINDOW_POP_PUSH × 16 = 32` 行，实现窗口推进：
$$
\text{win\_start}_{k+1} = \text{win\_start}_k + \text{pop\_push\_rows}
$$

## 3. 从Tile到解码器输入

对于 tile 内每一行，代码需要对 256 维数据（128 历史 + 111 新字 + 16 parity + 1 overall）建立一一对应：

- **旧历史位（列 0..127）**：对每个 \(k \in [0, N-1]\)，根据当前行的全局坐标 \((R, r)\)（整 tile 行索引和子块内行索引），先求
  $$
  br = (R \oplus 1) - 2G - 2\frac{N}{B} + 2\left\lfloor \frac{k}{B} \right\rfloor,\quad
  bc = \left\lfloor \frac{k}{B} \right\rfloor
  $$
  并得到块内偏移
  $$
  rr_{\text{blk}} = (k \bmod B) \oplus r,\quad
  cc_{\text{blk}} = r
  $$
  最终映射到矩阵坐标
  $$
  rr = br \cdot B + rr_{\text{blk}},\quad
  cc = bc \cdot B + cc_{\text{blk}},
  $$
  将 `Lch + La` 填入到 `lin_matrix[row_idx][k]`。

- **新系统位（列 128..238）及 parity/overall（列 239..255）**：位于当前 tile 的局部矩阵，可直接访问。以系统位 \(k = N + i\)（\(i \in [0, 127]\)）为例，
  $$
  C_t = \left\lfloor \frac{k - N}{B} \right\rfloor,\quad
  c_t = ((k - N) \bmod B) \oplus r
  $$
  对应列索引
  $$
  col = C_t \cdot B + c_t,
  $$
  parity（\(k = K + j\)）和 overall（\(k = BCH\_OVERALL\_IDX\)）使用同样的公式。由于这些列就在当前 tile 内，所以直接读取并写入 `lin_matrix[row_idx][k]`。

最终得到输入解码器的LLR矩阵： `lin_matrix`：包含信道 LLR 与外信息之和 \(L_\text{in} = L_\text{ch} + L_\text{a}\)；

## 4. 码芯处理与早停

- `tile_should_early_stop()` 会把 `lin_matrix` 做硬判决并运行一次 BCH(255,239) + overall parity，若全部满足则记为“tile 早停”，统计在 `TileEarlyStopCounter`。
- Tile 会调用 `Decoder_Core_(plain)`来进行BCH Chase 解码，核心步骤为：
  1. 将每行 256 LLR 拷入 `LinVec` 和 `LchVec`；
  2. 若该Tile是硬判Tile，走硬判决回退（`perform_hard_decode()`）；
  3. 否则调用 `chase_decode_256_(plain)`，输出仅包含外信息的 \( \omega_j \)矩阵；

## 5. 外信息归一化与写回

在软判决模式下，`process_tile_impl()` 会在得到仅包含外信息的 \( \omega_j \)矩阵后执行两级缩放：

1. **归一化**：若开启 `normalize_extrinsic`，统计所有有效输出$\omega_j$ 的平均绝对值幅值 \(g_\alpha\)，然后将 \( \omega_j \gets \omega_j / g_\alpha \)。
2. **全局 α**：再乘以 `tile_params.ALPHA`（每个Tile的alpha参数不同），对应 Pyndiah (21) 式：
   $$
   L_\text{next} = L_\text{ch} + \alpha \cdot \omega
   $$

写回阶段遵循与编码完全对偶的寻址公式，把 256 维外信息分别填回。

## 6. 输出组合

窗口全部处理完后，`ofec_decode_llr_impl()` 将 `channel_llr` 与 `work_llr` 相加：
$$
L_\text{post}(r,c) = L_\text{ch}(r,c) + L_\text{extr}(r,c)
$$
作为 `post_decoder_llr` 返回。`run_pipeline()` 随后用 `rx_info_from_bit_llr()` 抽取信息位，与发射端参考比较 BER。

通过以上步骤，oFEC 解码实现了“窗口滑动 + tile 迭代 + Chase SISO” 的完整流水线。
