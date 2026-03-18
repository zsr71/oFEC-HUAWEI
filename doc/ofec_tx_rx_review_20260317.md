# oFEC 发端/收端检查记录

日期：2026-03-17

这份记录总结了本次对 oFEC 发端、收端和 BER 统计链路的检查结果。重点关注“为什么当前纠后 BER 看起来偏低”。

## 结论概览

当前结果里至少有三类问题会直接影响你看到的 BER：

1. BER 参考序列 `tx_info_bits_ref` 不是原始 `info_bits` 的真实逆映射结果，但这一点更像“指标口径差异”，不一定是解码算法 bug。
2. BER 统计额外裁掉了大量前后比特，且打印出来的 `cut` 数量与真实裁剪量不一致。
3. 收端最终 `post_decoder_llr` 的组合方式与 `work_llr` 写回路径不完全一致，导致部分译码输出没有进入最终 BER 口径。

另有两个额外风险：

4. 已知零前缀只被赋值成较小的正 LLR，和当前早停阈值量级不一致，可能影响收敛和早停判断。
5. `detect_mode=2` 的 early-stop 只检查 LLR 可靠性，但命中后直接走 `row_early_stop_process_1()`，会继续给当前符号加 `±beta`，可能把错误符号“钉死”。
6. float 与 qfloat 两条路径在解码核心里使用的数值量纲不同，但当前 `beta/threshold` 等参数是共用的，这会导致“量化后反而比浮点好”的假象，或至少造成不公平比较。
7. `ebchPF` 版 Chase 把 overall parity 位也纳入“最不可靠位”集合，但 BCH 候选生成只基于前 255 位核心码字，这会浪费测试向量，甚至错误淘汰候选。
8. 在 `CHASE_SBR=2` 且 `SISO` 预算不足 32 的场景里，送入 Chase 的 32 行顺序与“低索引优先保留”语义不一致，可能把预算给到较旧的那半行。
9. hard tile 分支里 `perform_hard_decode()` 计算外信息时减掉的是 `Lin256` 而不是 `Lch256`，一旦启用 hard tile 会产生错误的 extrinsic。

## 问题 1：信息位提取口径不是源比特真逆映射

### 现象

`run_pipeline()` 在发端把编码矩阵转成 `tx_llr_mat` 后，调用 `rx_info_from_bit_llr()` 生成参考序列 `tx_info_bits_ref`，后续 pre/post BER 都拿它当真值：

- `src/common/pipeline/pipeline_runner.cpp:122-129`
- `src/common/pipeline/pipeline_runner.cpp:255-274`

但 `rx_info_from_bit_llr()` 当前实现只是：

1. 跳过前 `tile_height_rows()` 行；
2. 然后逐行顺序抽取前 `111` 列。

对应代码：

- `src/rx/extract/info_extract.cpp:31-49`

### 为什么有问题

发端编码时，原始 `info_bits` 并不是按“矩阵行内前 111 列顺序”直接写进去的。

真实映射分两层：

1. 先通过 `u_index(R, r, k)` 把原始比特流映射到每个编码行的 `111` 个系统位位置。
2. 再通过 `ct = (k % B) ^ r` 把系统位写回矩阵列坐标。

对应代码：

- `src/tx/ofec/ofec_encoder.cpp:60-75`
- `src/tx/ofec/ofec_encoder.cpp:77-94`
- `src/tx/ofec/ofec_encoder.cpp:110-141`

因此，`rx_info_from_bit_llr()` 现在并不是 `ofec_encode()` 的逆映射。

### 本地验证结果

我用当前仓库代码做了两个验证：

1. 直接比较 `info_bits` 与当前 `tx_info_bits_ref`
2. 用和 `ofec_encode()` 完全一致的 `u_index()` 逻辑反解编码矩阵

结果：

- 当前实现下，`info_bits` 与 `tx_info_bits_ref` 有 `3,107,958 / 6,251,520` 位不一致。
- 用编码器同构逆映射反解后，可做到 `0 / 6,251,520` mismatch。

这说明 `rx_info_from_bit_llr()` 不是“原始 `info_bits`”意义上的逆映射。

### 对当前目标的影响怎么理解

这里需要区分两个目标：

1. 如果目标是验证“解码后矩阵经过同一提取器后，是否能回到发端矩阵同一提取器下的结果”，那么当前口径是自洽的。
2. 如果目标是验证“原始发送 `info_bits` 的真实 BER”，那么当前口径不是源比特 BER。

也就是说：

- 高信噪比下 BER 仍然完全可能到 0；
- 这并不能证明 `rx_info_from_bit_llr()` 是源比特逆映射；
- 只能说明“发端参考 + 收端提取”这套闭环在当前口径下是一致的。

### 影响

- 对“算法是否把矩阵恢复正确”的验证，这一项不是首要 bug。
- 对“原始源比特 BER”的解释，这一项会造成口径偏差。
- 如果后续要和 baseline 做“源比特 BER”对齐，这里仍然需要修。

## 问题 2：BER 统计裁剪区间过大，且日志少报

### 现象

`compute_ber()` 当前不是全量比较，而是会跳过：

- 前 `4 * win_bits`
- 后 `1 * win_bits`

对应代码：

- `src/rx/ber/ber.cpp:18-27`

但注释写的是“去掉首尾各一个 window 覆盖的比特”，和实际实现不一致。

此外，`compute_and_print_ber()` 打印 `cut_total` 时只按“前 1 个 window + 后 1 个 window”计算，和真实裁剪量不一致：

- `src/rx/ber/ber.cpp:54-64`

### 本地验证结果

按默认参数：

- `win_bits = 156,288`
- `skip_prefix = 625,152`
- `skip_suffix = 156,288`
- 总裁剪量 = `781,440`
- 实际参与比较的比特数 = `5,470,080`

也就是说，当前 BER 统计直接不看约 `12.5%` 的数据，而且日志里的 `cut=` 还会少报。

### 影响

- 会系统性压低 BER。
- 不同参数组合下，统计口径不可直观看清。
- 你很容易误以为译码效果比实际更好。

## 问题 3：`post_decoder_llr` 合成路径与 `work_llr` 不一致

### 现象

窗口处理时，`process_window_impl()` 会把每个 tile 的输出写回 `work_llr`：

- `src/rx/ofec/detail/ofec_window_impl.ipp:112-130`

但 frame 级最后生成 `post_decoder_llr` 时，用的不是 `work_llr`，而是：

- `channel_llr + last_tile_history_llr`

对应代码：

- `src/rx/ofec/detail/ofec_decode_impl.ipp:65-70`
- `src/rx/ofec/detail/ofec_decode_impl.ipp:107-113`

同时，`last_tile_history_accum` 的更新只发生在历史位 `k < 128` 这条写回路径里：

- `src/rx/ofec/detail/ofec_tile_writeback.ipp:106-157`

而新系统位 / parity / overall 虽然已经写进了 `tile_out`：

- `src/rx/ofec/detail/ofec_tile_writeback.ipp:64-104`

却没有同步进入 `last_tile_history_accum`。

### 为什么有问题

最终 BER 提取使用的是 `decode_result.post_decoder_llr`：

- `src/common/pipeline/pipeline_runner.cpp:255-257`

如果 `post_decoder_llr` 的来源不是完整的 `channel_llr + work_llr`，而只是 `channel_llr + last_tile_history_llr`，那么：

- 某些已经写回到 `work_llr` 的更新不会进入最终 BER；
- `post-FEC` 统计口径和 tile 内部真实输出不一致。

另外，当前窗口实现里只有“最后一个 soft tile”或 hard tile 才会把 history 捕获到 `last_tile_history_accum`：

- `src/rx/ofec/detail/ofec_window_impl.ipp:60-64`
- `src/rx/ofec/detail/ofec_window_impl.ipp:88-100`

在当前常见配置“所有 tile 都是 soft decode”下，这意味着只有 `t == last_soft_tile_idx` 的那个 tile 会参与 `last_tile_history_llr` 的最终合成，其他 soft tile 的 history 写回不会进入 `post_decoder_llr`。

### 影响

- 纠后 BER 可能被进一步“美化”或扭曲。
- 调试时看 `work_llr` 和看最终 BER 可能对不上。

## 问题 4：已知零前缀的 LLR 幅度过小

### 现象

`apply_known_zero_prefix()` 对已知零前缀写死的是：

- `bit0_llr = 2.0f`

对应代码：

- `src/rx/llr/llr_known_prefix.cpp:15-25`

而你当前单次配置里，早停相关阈值已经是：

- `kEarlyStopV2LlrAbsThreshold = 26.0f`

对应代码：

- `apps/ofec_single.cpp:16-24`

### 为什么有问题

注释里说“强制为比特 0（大正 LLR）”，但 `+2.0` 并不算“大正”，至少和当前早停阈值量级完全不在一个范围里。

这不一定直接导致 BER 偏低，但会导致：

- warm-up 区域置信度不够强；
- 早停判断和实际“已知比特”假设不一致；
- 收敛行为和预期存在偏差。

### 影响

- 影响早停命中率和外信息传播稳定性。
- 会增加定位性能问题时的干扰项。

## 问题 5：`v2` early-stop 与 `process_1` 的组合有较高误停风险

### 现象

`detect_tile_early_stop_v2()` 的判据只看一行里“不可靠 bit”的个数：

- `src/rx/ofec/earlystop/tile_early_stop_stats2.ipp:22-59`

只要 `|LLR| < threshold` 的 bit 数不超过 `max_unreliable_bits`，该行就会被标成 `row_passed_flags[r] = true`。

但行级真正命中 early-stop 后，当前代码固定走的是 `row_early_stop_process_1()`：

- `src/rx/ofec/ofec_row_decoder_core.cpp:130-140`

而不是 `row_early_stop_process_2()`：

- `src/rx/ofec/earlystop/row_early_stop_process_2.ipp:10-20`

`process_1()` 的输出公式是：

- `Y2 = (Lin - Lch) + sign(Lin) * beta`

对应代码：

- `src/rx/ofec/earlystop/row_early_stop_process_1.ipp:16-20`

### 为什么有问题

`v2` 判据并不验证该行是否已经满足 BCH syndrome / overall parity，只验证“当前看起来够可靠”。

在这种前提下，一旦一行被标成 early-stop：

1. 不再跑 Chase；
2. 直接按 `sign(Lin)` 给外信息再加一层 `±beta`；
3. 会进一步强化当前符号方向。

如果这行其实只是“看起来很可靠，但还不是正确码字”，这条路径会把错误符号继续加固，而不是留给 Chase 去纠正。

### 影响

- 对 `detect_mode=2` 特别敏感；
- 在中高 SNR、少量强错误 bit 的情况下，可能出现“误停后锁错”的现象；
- 这类问题会直接把中段 BER 拉高，即使边界区统计已经被裁掉。

## 量化观察：当前 `clip_ratio=0.5` 属于激进配置

### 现象

动态 clip 的实现是：

- `clip = 第 ceil(ratio * N) 大的 |LLR|`

对应代码：

- `src/common/qfloat/qfloat.ipp:228-245`

文档里也明确写了，这样做会让大约 `ratio` 比例的样本进入饱和区：

- `doc/rx_llr_quantization.md:15-22`

### 说明

这不一定是代码 bug，但如果 `LLR_CLIP_RATIO = 0.5`，clip 会落在绝对值分布的中位附近，属于非常激进的量化配置。

在 `LLR_BITS = 6` 这类低位宽下，这种配置很容易：

- 压缩中高置信度 LLR 的动态范围；
- 让 Chase/early-stop 更难区分“可靠”和“非常可靠”；
- 直接损伤中段 BER。

因此它更像“高风险参数”，不是本次看到的明确实现 bug。

## 问题 6：float 与 quantized 路径使用了不同量纲，但参数没有分开标定

### 现象

qfloat 路径进入 Chase 核心时，并不是把“反量化后的真实 LLR 幅度”送进去，而是把 qfloat 的码值直接当作 `float` 使用：

- `src/rx/ofec/common/lin_matrix_adapters.ipp:20-33`

其中：

- `combine(qfloat, qfloat)` 返回 `Lch.code() + La.code()`
- `channel(qfloat)` 返回 `v.code()`

也就是说，量化路径的 `lin_matrix / lch_matrix` 在 core 里实际跑的是“码值域”。

而 float 路径进入 core 时，直接使用真实幅度：

- `src/rx/ofec/common/lin_matrix_adapters.ipp:5-18`

与此同时，Chase / early-stop 共用的参数 `beta`、`ALPHA`、`EARLY_STOP_V2_LLR_ABS_THRESHOLD` 都来自同一份 `Params`：

- `src/ofec_single/ofec_single_params.cpp:42-74`
- `apps/ofec_single.cpp:24-44`

### 为什么这会导致“量化后反而更好”

当前配置下，量化链路的动态 `clip` 是按每帧实时估计的；我用当前单次配置做了一个本地探针，得到：

- `clip ≈ 0.999972`
- `|LLR|` 中位数约 `0.999972`
- `|LLR| p99 ≈ 2.09142`
- 6 bit qfloat 的 `Q = 31`

在这个量纲下：

- quantized 路径里，码值 `31` 反量化后只对应约 `1.0`
- 码值 `26` 反量化后只对应约 `0.84`

但你当前参数里：

- `beta_list = {8.571428, 10.037715, 16.865997, 31.428572}`
- `EARLY_STOP_V2_LLR_ABS_THRESHOLD = 26.0`

这组数更像是“码值域参数”，而不是“float 幅度域参数”。

于是两条路径的实际含义变成：

1. 对 qfloat 路径：
   - `beta = 31` 大约对应“反量化后 1.0 左右”的外信息强度
   - `threshold = 26` 大约对应“反量化后 0.84 左右”的可靠度门限
2. 对 float 路径：
   - `beta = 31` 就是货真价实的 31 幅度外信息
   - `threshold = 26` 就是货真价实的 26 幅度门限

这两个量级差了一个数量级以上。

### 进一步放大差异的实现细节

qfloat 路径在每次 tile 输出后还会做一次“量化回目标精度”的硬裁剪：

- `src/rx/ofec/detail/ofec_tile_decode.ipp:102-114`
- `src/rx/ofec/common/lin_matrix_adapters.ipp:41-52`

而 float 路径这里完全不裁剪：

- `src/rx/ofec/common/lin_matrix_adapters.ipp:35-39`

因此：

- qfloat 路径天然带有“限幅/正则化”；
- float 路径则会把过大的 fallback / extrinsic 原样保留下来。

在当前 `beta_list` 很大的配置下，这会进一步放大“量化优于浮点”的现象。

### 影响

- float 与 quantized 的性能对比当前并不公平。
- 你看到的“量化后反而更好”很可能不是量化本身带来增益，而是：
  - qfloat 路径在码值域里运行；
  - 参数更接近码值域标定；
  - 同时还带有限幅保护。

### 建议

如果要严肃比较 float 和 quantized：

1. 先把 `beta`、`EARLY_STOP_V2_LLR_ABS_THRESHOLD` 这类参数区分为“float 幅度域”和“qfloat 码值域”两套；
2. 或者统一把 quantized 路径在 core 前先还原到 float 幅度域，再共用同一套参数；
3. 至少在当前代码结构下，不要直接拿同一组 `beta/threshold` 做 float vs qfloat 结论。

## 问题 7：`ebchPF` Chase 把 overall parity 位错误纳入不可靠位搜索

### 现象

`plain` 版 Chase 在选择最不可靠位时，只看前 `255` 个 core bit：

- `src/rx/ofec/plain/detail/chase256_plain_reliability.ipp:6-24`

但 `ebchPF` 版实现却把 `0..255` 共 `256` 位都纳入了最不可靠位集合，其中包含第 `255` 位 overall parity：

- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp:66-84`

后续生成测试向量时，也会对这些位置执行翻转：

- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp:263-266`

但真正的 BCH 硬译码调用只吃前 `255` 位：

- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp:268-272`
- `include/newcode/common/bch/bch_255_239.hpp:23-34`

### 为什么有问题

overall parity 位并不参与 BCH(255,239) 的核心候选搜索。

因此如果 `lrp_pos` 里选中了 `255` 号位：

1. 一部分 Chase pattern 会浪费在“翻转 overall parity”上；
2. 这些翻转对 255 位 BCH 候选本身没有搜索价值；
3. 等效上会减少真正用于 core bit 的搜索预算。

更糟的是，`ebchPF` 版后面还有一条针对 parity 的附加判据：

- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp:277-280`

如果 pattern 恰好翻了 `tmp_in[255]`，这里还有机会把本来有效的 BCH 候选额外判成 `ok = false`。

### 影响

- 这个问题会直接降低 `ebchPF` 版 Chase 的有效搜索能力。
- 在 `CHASE_L` 不大、`NTEST` 有限时，损失尤其明显。
- 这类问题会表现成：
  - 同样参数下 `ebchPF` 比预期更差；
  - 或者 `plain` / `ebchPF` 的性能差异异常。

## 问题 8：`CHASE_SBR=2` 时 row 索引顺序与低索引优先预算不一致

### 现象

`prepare_tile_inputs()` 在 `CHASE_SBR=2` 时，会把 tile 底部 32 行送入 row core / Chase。

当前构造顺序是：

- `row_idx = 0..15` 对应局部行 `320..335`
- `row_idx = 16..31` 对应局部行 `336..351`

对应代码：

- `src/rx/ofec/detail/ofec_tile_input.ipp:143-152`

也就是说，进入 Chase 的 32 行是“上半 16 行在前，下半 16 行在后”。

但 MUX / budget 裁剪在 `G=1` 路径下采用的是“低索引优先保留”：

- `src/rx/ofec/mux/mux_siso_budget.cpp:22-39`

因此当 `siso_active_for_tile < 32` 时，被优先保留的是 `row_idx` 较小的那一半，也就是 `320..335`，而不是更靠近 tile 底部的 `336..351`。

### 为什么有问题

如果你的设计意图是“优先保证最新、最靠下的 16 行先拿到 Chase 预算”，当前顺序与这个目标相反。

这不会影响：

- `siso_active_for_tile = 32`
- 或者 `CHASE_SBR = 1`

但会影响：

- `CHASE_SBR = 2`
- 且 `siso_active_for_tile < 32`
- 并且预算裁剪按低索引优先时

### 影响

- 对你当前 `SISO_ACTIVE_LIST = {32,32,32,32}` 配置，这一项不是首要问题。
- 但在后续做降 SISO、MUX 受限、grouped budget、reconfig 调度实验时，这会直接影响“哪 16 行先被 Chase”。
- 这类顺序偏差会让中段 BER 变差，而且很难只靠 aggregate 统计看出来。

## 问题 9：hard tile 分支的 extrinsic 计算用错了输入

### 现象

`perform_hard_decode()` 的接口同时接收了：

- `Lin256`
- `Lch256`

对应声明：

- `include/newcode/ofec_decoder_hard.hpp:9-14`

但实现里最后构造外信息时，减掉的并不是 `Lch256`，而是 `Lin256`：

- `src/rx/ofec/ofec_hard_decode_bch.cpp:40-45`

当前代码：

```cpp
const float Lpost = sign * hard_mag;
const float Lch = qfloat::llr_to_float(Lin256[i]);
Y2[i] = Lpost - Lch;
```

### 为什么有问题

hard decode 输出的应该是 extrinsic，语义上应当满足：

```text
extrinsic = Lpost - Lch(channel)
```

而不是：

```text
extrinsic = Lpost - Lin(total)
```

如果 `Lin = Lch + La`，那么当前实现会把已有的 a priori 信息也一起减掉，导致 hard 分支输出的外信息幅度和方向都不再是预期的 extrinsic。

### 影响

- 对你当前默认的全软配置，这一项通常不生效。
- 一旦后续启用 hard tile、混合软硬译码，结果会直接受影响。
- 它会让 hard decode 分支和 soft decode 分支的输出语义不一致，增加混合调度时的误码风险。

## 建议的修复顺序

建议按下面顺序处理，否则后续 BER 对比没有意义：

1. 先修 `compute_ber()` 的裁剪口径，并让日志打印真实裁剪量。
2. 再统一 `post_decoder_llr` 的组合方式，保证它和 `work_llr` 写回链路一致。
3. 然后重点排查 `detect_mode=2 + row_early_stop_process_1()` 这条 early-stop 路径是否过于激进。
4. 再修 `ebchPF` Chase 的不可靠位集合，确保 overall parity 不进入 `lrp_pos`。
5. 然后把 float/qfloat 的参数量纲统一，至少分开标定 `beta` 与 early-stop threshold。
6. 最后再调量化、alpha/beta、阈值等性能参数。
7. 如果后续要和 baseline 对齐“源比特 BER”，再单独修 `rx_info_from_bit_llr()` 的逆映射口径。
8. 如果后续要做受限 SISO / MUX 实验，再检查 `CHASE_SBR=2` 时 row 顺序是否需要改成“更靠底部的行优先”。
9. 如果后续要启用 hard tile，再修 `perform_hard_decode()` 里的 `Lch256` 使用错误。

## 相关代码位置总表

- `src/common/pipeline/pipeline_runner.cpp:122-129`
- `src/common/pipeline/pipeline_runner.cpp:255-274`
- `src/rx/extract/info_extract.cpp:31-49`
- `src/tx/ofec/ofec_encoder.cpp:60-75`
- `src/tx/ofec/ofec_encoder.cpp:77-94`
- `src/tx/ofec/ofec_encoder.cpp:110-141`
- `src/rx/ber/ber.cpp:18-27`
- `src/rx/ber/ber.cpp:54-64`
- `src/rx/ofec/detail/ofec_window_impl.ipp:112-130`
- `src/rx/ofec/detail/ofec_window_impl.ipp:60-64`
- `src/rx/ofec/detail/ofec_window_impl.ipp:88-100`
- `src/rx/ofec/detail/ofec_decode_impl.ipp:65-70`
- `src/rx/ofec/detail/ofec_decode_impl.ipp:107-113`
- `src/rx/ofec/detail/ofec_tile_writeback.ipp:64-104`
- `src/rx/ofec/detail/ofec_tile_writeback.ipp:106-157`
- `src/rx/llr/llr_known_prefix.cpp:15-25`
- `src/rx/ofec/earlystop/tile_early_stop_stats2.ipp:22-59`
- `src/rx/ofec/ofec_row_decoder_core.cpp:130-140`
- `src/rx/ofec/earlystop/row_early_stop_process_1.ipp:16-20`
- `src/rx/ofec/earlystop/row_early_stop_process_2.ipp:10-20`
- `src/common/qfloat/qfloat.ipp:228-245`
- `doc/rx_llr_quantization.md:15-22`
- `src/rx/ofec/common/lin_matrix_adapters.ipp:20-33`
- `src/rx/ofec/detail/ofec_tile_decode.ipp:102-114`
- `src/rx/ofec/common/lin_matrix_adapters.ipp:41-52`
- `src/rx/ofec/plain/detail/chase256_plain_reliability.ipp:6-24`
- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp:66-84`
- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp:263-280`
- `include/newcode/common/bch/bch_255_239.hpp:23-34`
- `src/rx/ofec/detail/ofec_tile_input.ipp:143-152`
- `src/rx/ofec/mux/mux_siso_budget.cpp:22-39`
- `include/newcode/ofec_decoder_hard.hpp:9-14`
- `src/rx/ofec/ofec_hard_decode_bch.cpp:40-45`
- `apps/ofec_single.cpp:16-24`
