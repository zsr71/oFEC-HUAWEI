# `new_float_only` 比特流转代码导读

本文按“发端比特 -> 编码 -> 调制 -> 信道 -> 解调 -> 解码 -> BER”的顺序，列出 [`new_float_only`](/home/zsr71/projects/newcode/new_float_only) 里最关键的代码文件，方便从整体上检查一条比特是如何走完整个链路的。

## 1. 从哪里开始看

如果你只想先抓主线，建议从这几个入口开始：

1. [`apps/ofec_single_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_single_float.cpp)
2. [`src/ofec_single/ofec_single_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/ofec_single/ofec_single_runner.cpp)
3. [`src/pipeline/pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

其中真正把整条链串起来的是：

- [`run_pipeline()`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

## 2. 总体数据流

当前 `new_float_only` 的主链路可以概括成：

```text
info_bits
-> oFEC encode
-> coded_bits
-> QAM/BPSK modulate
-> AWGN channel
-> soft demod (LLR)
-> LLR matrix
-> known prefix processing
-> plain oFEC decoder
-> post_decoder_llr
-> info extraction
-> BER
```

对应主线代码都在：

- [`src/pipeline/pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

## 3. App 层入口

### 3.1 单次运行入口

- [`apps/ofec_single_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_single_float.cpp)

作用：

- 构造 `SingleRunConfig`
- 设置 Eb/N0、seed、解码参数、trace 开关
- 调用 `run_single()`

### 3.2 单次运行包装层

- [`src/ofec_single/ofec_single_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/ofec_single/ofec_single_runner.cpp)

作用：

- 规范化 decoder 配置
- 创建日志
- 构造 `PipelineConfig`
- 调用 `run_pipeline()`

如果只关心“真正的链路”，可以从这里继续跳到：

- [`src/pipeline/pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

## 4. 发端：原始比特生成

关键文件：

- [`src/tx/bitgen/bitgen.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/bitgen/bitgen.cpp)

关键函数：

- `bitgen::generate_bits(num_bits, seed, random_bits)`

作用：

- 生成发送端原始信息比特 `info_bits`
- `random_bits=false` 时会全 0
- `random_bits=true` 时按 seed 生成伪随机比特

在主线里的调用位置：

- [`pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

## 5. 发端：oFEC 编码

关键文件：

- [`src/tx/ofec/ofec_encoder.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/ofec/ofec_encoder.cpp)

关键函数：

- `ofecencoder::ofec_encode(bits, params)`

作用：

- 把输入信息比特编码成二维 `code_matrix`
- 内部会构造历史部分、系统位部分、BCH parity 和 overall parity
- 输出是一个二维比特矩阵，而不是一维码流

这一步是当前发送端最核心的编码逻辑。

## 6. 发端：把编码矩阵展开成发送比特

关键代码位置：

- [`pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

相关语句：

- `flatten_row_major(code_matrix, ...)`

作用：

- 按行优先把 `code_matrix` 展平为一维 `coded_bits`
- 当前工程默认是“无交织直通”，所以这里没有额外 permutation

## 7. 发端：调制

关键文件：

- [`src/tx/qam/qam.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/qam/qam.cpp)

关键函数：

- `mod::qam_modulate(bits, n_bps)`

作用：

- `n_bps=1` 时走 BPSK
- `n_bps>1` 且为偶数时走 Gray QAM
- 输出复数符号序列 `tx_syms`

这一步把一维编码比特变成复平面发射符号。

## 8. 信道：加 AWGN

关键文件：

- [`src/channel/awgn.cpp`](/home/zsr71/projects/newcode/new_float_only/src/channel/awgn.cpp)
- [`src/channel/ebn0_to_sigma.cpp`](/home/zsr71/projects/newcode/new_float_only/src/channel/ebn0_to_sigma.cpp)

关键函数：

- `channel::add_awgn(x, ebn0_dB, bits_per_symbol, seed)`

作用：

- 根据 `Eb/N0` 和调制阶数计算噪声标准差
- 用 `channel_seed` 生成高斯噪声
- 输出加噪后的 `rx_syms`

## 9. 接收端：软解调得到 LLR

关键文件：

- [`src/rx/qam/qam_llr.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/qam/qam_llr.cpp)

关键函数：

- `demod::qam_llr_from_ebn0(y, n_bps, ebn0_dB, code_rate)`
- 内部调用 `qam_llr_logsumexp(...)`

作用：

- 把接收符号 `rx_syms` 转成一维 `llr`
- BPSK 直接取实部
- QAM 使用 log-sum-exp 精确 LLR

## 10. 接收端：LLR 回填成矩阵

关键文件：

- [`src/rx/llr/ofec_llr_matrix.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/llr/ofec_llr_matrix.cpp)

关键函数：

- `llr_to_matrix_row_major(llr, rows, cols)`

作用：

- 把一维 soft LLR 按行优先回填为 `llr_mat`
- 这个矩阵形状和发端 `code_matrix` 对齐

## 11. 接收端：已知前缀处理

关键文件：

- [`src/rx/llr/llr_known_prefix.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/llr/llr_known_prefix.cpp)

关键函数：

- `apply_known_zero_prefix(llr_mat, params)`

作用：

- 把已知前缀区域强制设成 bit0 对应的大正 LLR
- 如果 `NORMALIZE_KNOWN_PREFIX_TAIL=true`，还会对尾部区域做一次平均幅度归一化

这一步发生在真正解码之前。

## 12. 解码入口：plain float decoder

关键文件：

- [`src/rx/ofec/ofec_frame_decode.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/ofec_frame_decode.cpp)

关键函数：

- `decode_plain_llr(channel_llr, params, normalize_extrinsic, tx_llr_ref)`

作用：

- 这是整个 oFEC 解码器的顶层入口
- 真正的主体逻辑继续进入 detail 层

## 13. 解码主框架：window 级

关键文件：

- [`src/rx/ofec/detail/ofec_decode_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_decode_impl.ipp)

关键函数：

- `decode_plain_llr_impl(...)`

作用：

- 创建 `work_llr`
- 创建 `last_tile_history_llr`
- 按窗口滑动调用 `process_window_impl(...)`
- 最终输出：

```text
post_decoder_llr = channel_llr + last_tile_history_llr
```

这是当前解码框架最重要的总控文件。

## 14. 解码主框架：window 内 tile 调度

关键文件：

- [`src/rx/ofec/detail/ofec_window_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_window_impl.ipp)

关键函数：

- `process_window_impl(...)`

作用：

- 在一个 window 内逐个 tile 处理
- 从 `work_llr` / `channel_llr` 中切出当前 `tile_in` / `ch_tile`
- 为每个 tile 选择：
  - `ALPHA`
  - `beta`
  - 是否 hard decode
- 调用 `process_tile_impl(...)`
- 把 tile 输出写回 `work_llr`

## 15. 解码主框架：单个 tile 流程

关键文件：

- [`src/rx/ofec/detail/ofec_tile_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_impl.ipp)

关键函数：

- `process_tile_impl(...)`

作用：

1. `prepare_tile_inputs(...)`
2. `decode_tile(...)`
3. `writeback_tile(...)`

也就是 tile 级完整解码闭环。

## 16. Tile 输入重排：128 列 -> 256 维 Chase 输入

关键文件：

- [`src/rx/ofec/detail/ofec_tile_input.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_input.ipp)

关键函数：

- `prepare_tile_inputs(...)`

作用：

- 从 `tile_in` 和 `ch_tile` 里取出当前 tile 的数据
- 生成：
  - `lin_matrix`
  - `lch_matrix`
- 其中：
  - `k < 128` 部分对应历史区
  - `k >= 128` 部分对应当前行的系统位、校验位和 overall parity

这一步是“把 oFEC tile 布局映射到 Chase 256 维行码字”的核心。

## 17. Tile 内部译码：Chase / hard decode

关键文件：

- [`src/rx/ofec/detail/ofec_tile_decode.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_decode.ipp)
- [`src/rx/ofec/plain/chase256_plain.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/plain/chase256_plain.cpp)
- [`src/rx/ofec/ofec_hard_decode_bch.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/ofec_hard_decode_bch.cpp)

关键过程：

1. `decode_tile(...)` 调用 row core
2. soft 路径进入 `chase::Decoder_Core_plain`
3. hard 路径进入 `perform_hard_decode(...)`

作用：

- soft 路径输出 256 维 `extrinsic`
- 可选做 `normalize_extrinsic`
- 再乘当前 tile 的 `ALPHA`

如果你要进一步看 Chase 算法本体，可以继续追：

- [`src/rx/ofec/plain/detail/chase256_plain_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/plain/detail/chase256_plain_impl.ipp)

## 18. Tile 输出写回

关键文件：

- [`src/rx/ofec/detail/ofec_tile_writeback.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_writeback.ipp)

关键函数：

- `writeback_tile(...)`

作用：

- 把 256 维行级 `extrinsic` 写回到 tile 的 128 列布局
- 更新当前 `tile_out`
- 必要时累计 `last_tile_history_accum`

这一步把 Chase 输出重新映射回全局工作矩阵。

## 19. 后处理：从后验 LLR 中提取信息位

关键文件：

- [`src/rx/extract/info_extract.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/extract/info_extract.cpp)

关键函数：

- `matrix::rx_info_from_bit_llr(bit_llr_mat, params)`

作用：

- 从最终 LLR 矩阵里提取信息位区域
- 做硬判决
- 返回一维 `rx_info_bits`

这一步在主线里会做两次：

- 一次针对 `llr_mat`，得到 pre-FEC
- 一次针对 `post_decoder_llr`，得到 post-FEC

## 20. BER 统计

关键文件：

- [`src/rx/ber/ber.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ber/ber.cpp)

关键函数：

- `compute_ber(...)`
- `compute_and_print_ber(...)`

作用：

- 比较发端参考信息位和接收端硬判信息位
- 当前实现会裁掉前后若干个 window 覆盖区域
- 输出：
  - `errors`
  - `total`
  - `ber`

## 21. 如果只想顺着一条主线读

最推荐的阅读顺序是：

1. [`apps/ofec_single_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_single_float.cpp)
2. [`src/ofec_single/ofec_single_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/ofec_single/ofec_single_runner.cpp)
3. [`src/pipeline/pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)
4. [`src/tx/bitgen/bitgen.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/bitgen/bitgen.cpp)
5. [`src/tx/ofec/ofec_encoder.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/ofec/ofec_encoder.cpp)
6. [`src/tx/qam/qam.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/qam/qam.cpp)
7. [`src/channel/awgn.cpp`](/home/zsr71/projects/newcode/new_float_only/src/channel/awgn.cpp)
8. [`src/rx/qam/qam_llr.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/qam/qam_llr.cpp)
9. [`src/rx/llr/ofec_llr_matrix.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/llr/ofec_llr_matrix.cpp)
10. [`src/rx/llr/llr_known_prefix.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/llr/llr_known_prefix.cpp)
11. [`src/rx/ofec/ofec_frame_decode.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/ofec_frame_decode.cpp)
12. [`src/rx/ofec/detail/ofec_decode_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_decode_impl.ipp)
13. [`src/rx/ofec/detail/ofec_window_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_window_impl.ipp)
14. [`src/rx/ofec/detail/ofec_tile_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_impl.ipp)
15. [`src/rx/ofec/detail/ofec_tile_input.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_input.ipp)
16. [`src/rx/ofec/detail/ofec_tile_decode.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_decode.ipp)
17. [`src/rx/ofec/detail/ofec_tile_writeback.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_writeback.ipp)
18. [`src/rx/extract/info_extract.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/extract/info_extract.cpp)
19. [`src/rx/ber/ber.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ber/ber.cpp)

## 22. 一句话总结

如果只抓一句主线，可以记成：

> `run_pipeline()` 把发送比特生成、oFEC 编码、调制、加噪、软解调和 plain oFEC 解码串成一条链；解码内部再由 `decode_plain_llr_impl()` 组织成 `window -> tile -> Chase core -> writeback -> posterior LLR` 的流程。
