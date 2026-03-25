# `new_float_only` 数据对象流转图

本文基于 [`BIT_FLOW_CODE_GUIDE.md`](/home/zsr71/projects/newcode/new_float_only/BIT_FLOW_CODE_GUIDE.md)，专门从“数据对象”的角度梳理一次完整链路里每个关键变量是怎么流动的。

## 1. 总览图

```text
info_bits
  std::vector<uint8_t>
  发送端原始信息比特
        |
        v
code_matrix
  matrix::Matrix<uint8_t>
  oFEC 编码后的二维比特矩阵
        |
        v
tx_llr_mat
  matrix::Matrix<float>
  由 code_matrix 映射出的发送端参考 LLR 矩阵
        |
        +----------------------------+
        |                            |
        |                            v
        |                     tx_info_bits_ref
        |                     std::vector<uint8_t>
        |                     用于最终 BER 对比的发送端参考信息位
        |
        v
coded_bits
  std::vector<uint8_t>
  按行优先展开后的发射比特流
        |
        v
tx_syms
  std::vector<std::complex<float>>
  调制后的复数符号
        |
        v
rx_syms
  std::vector<std::complex<float>>
  加 AWGN 后的接收符号
        |
        v
llr
  std::vector<float>
  软解调得到的一维 LLR
        |
        v
llr_mat
  matrix::Matrix<float>
  按行优先回填后的接收 LLR 矩阵
        |
        +----------------------------+
        |                            |
        |                            v
        |                     rx_info_bits_pre
        |                     std::vector<uint8_t>
        |                     解码前硬判信息位
        |
        v
post_decoder_llr
  matrix::Matrix<float>
  plain oFEC 解码后的后验 LLR
        |
        v
rx_info_bits_post
  std::vector<uint8_t>
  解码后硬判信息位
        |
        +----------------------------+
        |                            |
        v                            v
pre_fec BER                    post_fec BER
```

## 2. 主链路变量流

最核心的主链路都在：

- [`pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

按代码顺序看，关键对象依次是：

### 2.1 `info_bits`

类型：

- `std::vector<uint8_t>`

来源：

- `bitgen::generate_bits(...)`

代码位置：

- [`bitgen.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/bitgen/bitgen.cpp)
- [`pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

含义：

- 最原始的发送信息比特

### 2.2 `code_matrix`

类型：

- `matrix::Matrix<uint8_t>`

来源：

- `ofecencoder::ofec_encode(info_bits, params)`

代码位置：

- [`ofec_encoder.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/ofec/ofec_encoder.cpp)

含义：

- oFEC 编码后的二维矩阵
- 后面调制前会先按行优先展开

### 2.3 `tx_llr_mat`

类型：

- `matrix::Matrix<float>`

来源：

- `hard_bits_to_llr_matrix(code_matrix, tx_ref_llr)`

代码位置：

- [`pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

含义：

- 把发端编码比特映射成一个理想参考 LLR 矩阵
- 主要用来提取发送端参考信息位，以及可选 Chase trace 对照

### 2.4 `tx_info_bits_ref`

类型：

- `std::vector<uint8_t>`

来源：

- `matrix::rx_info_from_bit_llr(tx_llr_mat, params)`

代码位置：

- [`info_extract.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/extract/info_extract.cpp)

含义：

- 从发送端参考矩阵中抽取的信息位
- 这是后续 BER 比较时的“真值”

### 2.5 `coded_bits`

类型：

- `std::vector<uint8_t>`

来源：

- `flatten_row_major(code_matrix, ...)`

代码位置：

- [`pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)

含义：

- 发往调制器的一维码流
- 当前实现没有额外交织

### 2.6 `tx_syms`

类型：

- `std::vector<std::complex<float>>`

来源：

- `mod::qam_modulate(coded_bits, n_bps)`

代码位置：

- [`qam.cpp`](/home/zsr71/projects/newcode/new_float_only/src/tx/qam/qam.cpp)

含义：

- 复数调制符号序列

### 2.7 `rx_syms`

类型：

- `std::vector<std::complex<float>>`

来源：

- `channel::add_awgn(tx_syms, ebn0_dB, n_bps, channel_seed)`

代码位置：

- [`awgn.cpp`](/home/zsr71/projects/newcode/new_float_only/src/channel/awgn.cpp)

含义：

- 通过 AWGN 信道后的接收符号

### 2.8 `llr`

类型：

- `std::vector<float>`

来源：

- `demod::qam_llr_from_ebn0(rx_syms, n_bps, ebn0_dB, code_rate)`

代码位置：

- [`qam_llr.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/qam/qam_llr.cpp)

含义：

- 接收端软解调输出的一维 LLR

### 2.9 `llr_mat`

类型：

- `matrix::Matrix<float>`

来源：

- `llr_to_matrix_row_major(llr, code_matrix.rows(), code_matrix.cols())`

代码位置：

- [`ofec_llr_matrix.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/llr/ofec_llr_matrix.cpp)

含义：

- 把一维 LLR 回填到二维矩阵
- 作为解码器输入

### 2.10 `rx_info_bits_pre`

类型：

- `std::vector<uint8_t>`

来源：

- `matrix::rx_info_from_bit_llr(llr_mat, params)`

含义：

- 解码前，直接从接收 LLR 上硬判得到的信息位

### 2.11 `post_decoder_llr`

类型：

- `matrix::Matrix<float>`

来源：

- `decode_plain_llr(llr_mat, params, normalize_extrinsic, &tx_llr_mat)`

代码位置：

- [`ofec_frame_decode.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/ofec_frame_decode.cpp)
- [`ofec_decode_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_decode_impl.ipp)

含义：

- plain oFEC 解码器输出的后验 LLR 矩阵

### 2.12 `rx_info_bits_post`

类型：

- `std::vector<uint8_t>`

来源：

- `matrix::rx_info_from_bit_llr(post_decoder_llr, params)`

含义：

- 解码后硬判得到的信息位

### 2.13 `pre_fec` / `post_fec`

类型：

- `BerStats`

来源：

- `compute_and_print_ber(tx_info_bits_ref, rx_info_bits_pre/post, ...)`

代码位置：

- [`ber.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ber/ber.cpp)

含义：

- 解码前后的 BER 统计结果

## 3. 解码器内部对象流

解码器内部的主线是：

```text
llr_mat
-> channel_llr
-> work_llr
-> tile_in / ch_tile
-> lin_matrix / lch_matrix
-> decoder_res.lout
-> tile_out
-> last_tile_history_llr
-> post_decoder_llr
```

## 4. 解码器内部总览图

```text
llr_mat
  matrix::Matrix<float>
  解码器输入
        |
        v
channel_llr
  matrix::Matrix<float>
  输入副本，作为固定信道项
        |
        +----------------------------------+
        |                                  |
        v                                  |
work_llr                                   |
  matrix::Matrix<float>                    |
  窗口内不断被 tile 输出覆盖               |
        |                                  |
        v                                  |
tile_in / ch_tile                          |
  matrix::Matrix<float>                    |
  当前 tile 的 a priori / channel 切片     |
        |                                  |
        v                                  |
lin_matrix / lch_matrix                    |
  matrix::Matrix<float>                    |
  送入 Chase core 的 256 维行输入          |
        |                                  |
        v                                  |
decoder_res.lout                           |
  matrix::Matrix<float>                    |
  当前 tile 的 256 维 extrinsic            |
        |                                  |
        v                                  |
tile_out                                   |
  matrix::Matrix<float>                    |
  写回当前 tile 布局后的输出               |
        |                                  |
        +----------------------+           |
        |                      |           |
        v                      v           |
work_llr updated      last_tile_history_llr|
                             累积最终 history
                                     |
                                     v
post_decoder_llr = channel_llr + last_tile_history_llr
```

## 5. `channel_llr`

创建位置：

- [`ofec_decode_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_decode_impl.ipp)

含义：

- `llr_mat` 的一个副本
- 表示固定不变的信道项
- 最终输出时会和 history 相加

## 6. `work_llr`

创建位置：

- [`ofec_decode_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_decode_impl.ipp)

含义：

- 窗口级工作矩阵
- 初始全 0
- 每处理完一个 tile，就把 tile 输出覆盖回去

它不是最终 posterior，而是迭代中的“工作状态”

## 7. `last_tile_history_llr`

创建位置：

- [`ofec_decode_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_decode_impl.ipp)

含义：

- 单独保存最终用于输出的 history
- 最终结果是：

```text
post_decoder_llr = channel_llr + last_tile_history_llr
```

## 8. `tile_in` / `ch_tile`

创建位置：

- [`ofec_window_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_window_impl.ipp)

含义：

- `tile_in`：当前 tile 的工作输入
- `ch_tile`：当前 tile 的信道输入

来源：

- 从 `work_llr` 和 `channel_llr` 中切片得到

## 9. `lin_matrix` / `lch_matrix`

创建位置：

- [`ofec_tile_input.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_input.ipp)

含义：

- `lin_matrix`：送入 Chase core 的总输入
- `lch_matrix`：只保留信道项的输入

关系：

```text
lin_matrix = Lch + La
lch_matrix = Lch
```

这是 tile 布局到 256 维行码字布局的关键重排结果。

## 10. `decoder_res.lout`

来源：

- `decode_tile(...)`
- 最终来自 `chase::Decoder_Core_plain(...)` 或 `perform_hard_decode(...)`

代码位置：

- [`ofec_tile_decode.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_decode.ipp)

含义：

- 当前 tile 内部 decoder 输出的 256 维 `extrinsic`

后处理：

1. 可选 `normalize_extrinsic`
2. 乘 `ALPHA`

## 11. `tile_out`

来源：

- `writeback_tile(...)`

代码位置：

- [`ofec_tile_writeback.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_writeback.ipp)

含义：

- 把 `decoder_res.lout` 重新映射回当前 tile 的 128 列布局之后的结果

后续用途：

- 覆盖回 `work_llr`

## 12. 一个最值得盯住的函数链

如果你要手动跟一条 bit 在解码器内部是怎么走的，最值得看的调用链是：

1. [`decode_plain_llr_impl()`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_decode_impl.ipp)
2. [`process_window_impl()`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_window_impl.ipp)
3. [`process_tile_impl()`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_impl.ipp)
4. [`prepare_tile_inputs()`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_input.ipp)
5. [`decode_tile()`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_decode.ipp)
6. [`writeback_tile()`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_writeback.ipp)

## 13. 一句话总结

如果只记一句数据对象流转，可以记成：

> 发端把 `info_bits` 编成 `code_matrix`，调制成 `tx_syms`，经信道后变成 `llr_mat`；解码器再把 `llr_mat` 通过 `work_llr -> tile_in -> lin_matrix -> lout -> tile_out -> last_tile_history_llr` 这一串对象转换，最终得到 `post_decoder_llr` 并提取出 `rx_info_bits_post`。
