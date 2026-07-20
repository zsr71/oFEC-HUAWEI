# 接收端：LLR 量化流程

本说明梳理 `run_pipeline()` 里浮点 LLR 在送入解码器前如何按 `Params::LLR_BITS` 做线性量化，并概述 `qfloat` 的编码方式。
## 1. 模式选择

`pick_llr_mode()` 根据 `params.LLR_BITS` 决定：

- `LLR_BITS = 16` → 直接以 `float` 传给解码器；
- `2 ≤ LLR_BITS ≤ 15` → 构造 `qfloat<LLR_BITS>`，对应码值区间 \([-Q, +Q]\)。

后续所有量化参数都打包进 `DecodeRequest`。

## 2. 动态裁剪幅度

当 `LLR_BITS < 16` 时，需要确定线性量化区间 \([-clip, +clip]\)。`run_pipeline()` 提供两种策略：

1. 如果 `Params::LLR_CLIP_RATIO > 0`，先将当前帧 LLR 的绝对值集合 \(\{|L_i|\}\) 按降序排列，取排名  
   $$
   k = \left\lceil \text{ratio} \cdot N \right\rceil,\quad clip = |L|_{(k)}
   $$
   确保只有`ratio` 比例的样本会被饱和；
`clip` 最终由 `DecodeRequest::quant_clip` 传递给具体解码器。

## 3. 线性量化与码值映射

`decoder_factory` 在收到量化请求后调用：
$$
code = \operatorname{sat}\!\left( \operatorname{round}\left( x \cdot \frac{Q}{clip} \right) \right),
\quad Q = 2^{(b-1)} - 1
$$
其中 `sat` 表示裁剪到 \([-Q, +Q]\)，`b = LLR_BITS`。`qfloat<b>` 以 16 bit 有符号整型保存 `code`，并提供：
$$
x \approx code \cdot \frac{clip}{Q}
$$
用于解码结束后的反量化。

## 4. 解码端使用方式

以 `ebchPF` 解码器为例：

1. 将接收 LLR 量化到 `Matrix<qfloat<b>>`；
2. `result.pre_decoder_llr` = 反量化后的浮点矩阵，仅用于日志/比较；
3. `ofec_decode_llr_ebchPF()` 直接在 `qfloat` 矩阵上完成窗口迭代；
4. 输出 `result.post_decoder_llr` 时再做一次反量化（或简单 cast）。

由于 `qfloat` 在乘法/加法时自动回到码值域，并且 `ALPHA`、`beta` 等缩放都在浮点空间计算完再重新量化，可以保证 oFEC Chase 解码的数值路径与纯浮点调试一致，只是幅值被线性压缩。
