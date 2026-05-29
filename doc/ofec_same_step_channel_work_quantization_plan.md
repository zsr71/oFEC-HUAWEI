# oFEC：Channel / Work 同步长分离量化改造方案

## 1. 背景

当前代码里，接收端量化参数只有一套：

- `Params::LLR_BITS`
- `Params::LLR_CLIP`
- `Params::LLR_CLIP_RATIO`

它们同时作用于：

1. `channel_llr`
2. `work_llr`
3. tile 内部 `tile_in / tile_out`
4. Chase core 的输入与输出量化

这会带来一个问题：

- 当 `channel` 被量化到较低位宽并发生饱和时，`work` 也被同样的位宽和范围限制；
- 这样即使 `work` 想在后续迭代中“压过”一个强负的 `channel`，可表示范围也不够。

本方案的目标是：

- 保持 `work` 和 `channel` 的量化步长完全一样；
- 只让 `work` 的可表示范围更大；
- 这样后续做加法、比较、判符号时，仍然可以把两边的 `code` 理解为同一个物理单位。

## 2. 目标设计

### 2.1 设计目标

希望支持如下配置：

- `channel_bits = 6`
- `work_bits = 7`
- `channel_clip = 8`
- `work_clip` 自动推导，使 `work` 与 `channel` 的量化步长严格一致

### 2.2 量化步长

对位宽 `b`，记：

```text
Q(b) = 2^(b - 1) - 1
```

则当前量化步长定义为：

```text
delta = clip / Q
```

若要求 `channel` 与 `work` 的步长完全相同，则必须满足：

```text
channel_clip / Qc = work_clip / Qw
```

其中：

- `Qc = Q(channel_bits)`
- `Qw = Q(work_bits)`

因此：

```text
work_clip = channel_clip * Qw / Qc
```

### 2.3 例子

当：

- `channel_bits = 6`，`Qc = 31`
- `work_bits = 7`，`Qw = 63`
- `channel_clip = 8`

则：

```text
work_clip = 8 * 63 / 31 = 16.2580645...
```

对应关系：

- `channel` 真实范围：`[-8, +8]`
- `work` 真实范围：`[-16.2580645, +16.2580645]`
- 两者量化步长完全相同：

```text
delta = 8 / 31 = 16.2580645 / 63
```

这个设计的直接好处是：

1. `work` 的真实幅度范围变大；
2. `work` 的量化精度不变；
3. `channel.code()` 与 `work.code()` 代表的是同一个“真实 LLR 单位格子”。

## 3. 为什么这个设计适合当前代码

当前 `qfloat` 路径里，tile 输入重排时的合成逻辑在：

- `src/rx/ofec/common/lin_matrix_adapters.ipp`

对量化路径的特化是：

```cpp
return static_cast<float>(Lch.code() + La.code());
```

这段逻辑成立的前提是：

1. `channel` 和 `work` 的量化步长相同；
2. 两边 `code` 对应同一个真实单位。

因此，如果我们采用“同步长、扩范围”的设计，那么：

- `channel.code()` 和 `work.code()` 仍然可以直接相加；
- 不需要把主路径改成“每次都先反量化到 float 再相加”；
- 现有 code-domain 的 Chase 输入口径可以继续保留。

这就是本方案相比“完全拆成不同实数口径再相加”的最大优势：

- 算法语义清晰；
- 改动可控；
- 运算实现更贴近现有框架。

## 4. 目标行为

### 4.1 旧行为

当前只有一套量化口径。例如：

- `channel_bits = 6`
- `work_bits = 6`
- `channel_clip = work_clip = 8`

则：

- `channel.code` 范围 `[-31, +31]`
- `work.code` 范围 `[-31, +31]`
- `Lin = channel.code + work.code`
- `Lin` 最多只能从 `-31` 被拉到 `0`

### 4.2 新行为

采用本方案后，例如：

- `channel_bits = 6`
- `work_bits = 7`
- `channel_clip = 8`
- `work_clip = 16.2580645`

则：

- `channel.code` 范围 `[-31, +31]`
- `work.code` 范围 `[-63, +63]`
- 两边步长相同
- `Lin = channel.code + work.code` 仍有物理意义

例如：

- `channel.code = -31`
- `work.code = +40`

则：

```text
Lin_code = +9
Lin_real = 9 * delta
```

这就允许 `work` 在后续迭代中真正压过一个已饱和的 `channel`。

## 5. 建议的参数改造

### 5.1 `ofec_single::Config`

当前：

- `llr_bits`
- `quant_clip_ratio`

建议拆成：

```cpp
std::size_t channel_llr_bits = 16;
std::size_t work_llr_bits = 16;
float channel_quant_clip_ratio = 0.0f;
bool work_keep_same_step_as_channel = true;
```

第一阶段不建议让用户直接手填 `work_clip`，而是：

1. 先根据 `channel_quant_clip_ratio` 算出 `channel_clip`
2. 再按 `work_keep_same_step_as_channel` 自动推导 `work_clip`

这样不容易把“同步长”关系配坏。

### 5.2 `Params`

当前：

```cpp
size_t LLR_BITS;
float  LLR_CLIP;
float  LLR_CLIP_RATIO;
```

建议拆成：

```cpp
size_t CHANNEL_LLR_BITS = 16;
float  CHANNEL_LLR_CLIP = 8.0f;
float  CHANNEL_LLR_CLIP_RATIO = 0.0f;

size_t WORK_LLR_BITS = 16;
float  WORK_LLR_CLIP = 8.0f;
bool   WORK_KEEP_SAME_STEP_AS_CHANNEL = true;
```

另外建议补两个便捷函数：

```cpp
static constexpr int qmax_from_bits(std::size_t bits);
float channel_llr_step() const;
float work_llr_step() const;
```

这样可以把步长关系集中在 `Params` 内部，不要散落在多个 `.cpp` 文件里重复计算。

### 5.3 `DecodeRequest`

当前：

```cpp
std::size_t quant_bits;
float quant_clip;
```

建议改为：

```cpp
std::size_t channel_quant_bits = 16;
float channel_quant_clip = 0.0f;
std::size_t work_quant_bits = 16;
float work_quant_clip = 0.0f;
```

理由：

1. `channel` 与 `work` 以后不再共享同一套位宽与范围；
2. decoder factory 需要同时知道两套量化口径；
3. 导出调试信息时，也需要明确区分“这是 channel 的 code 还是 work 的 code”。

## 6. 建议的数据流改造

### 6.1 当前数据流

当前量化路径大致是：

1. `channel_llr(float)` 量化成 `Matrix<qfloat<LLR_BITS>>`
2. `work_llr` 也是同一种 `qfloat`
3. `Lin = channel.code + work.code`
4. core 输出再按同一位宽写回 `work`

### 6.2 目标数据流

本方案建议改成：

1. `channel_llr(float)` 量化成 `Matrix<ChannelLLR>`
2. `work_llr` 单独存成 `Matrix<WorkLLR>`
3. `Lin = float(channel.code) + float(work.code)`，仍然在公共步长的 code-domain 中解释
4. core 输出 `lout` 继续按 code-domain `float` 表示
5. 写回 `work_llr` 时，按 `WorkLLR` 的范围截断和量化

这里的关键点是：

- `channel` 和 `work` 类型可以不同；
- 但只要步长相同，`code` 仍然可以直接相加。

## 7. 核心代码改造点

### 7.1 入口参数

涉及文件：

- `apps/ofec_single.cpp`
- `include/newcode/ofec_single_runner.hpp`
- `src/ofec_single/ofec_single_params.cpp`
- `src/ofec_single/ofec_single_summary.cpp`

需要做的事：

1. 把单一 `kLlrBits` 改成 `kChannelLlrBits` 和 `kWorkLlrBits`
2. 把单一 `kQuantClipRatio` 改成 `kChannelQuantClipRatio`
3. 在日志里打印：
   - `channel_bits`
   - `work_bits`
   - `channel_clip`
   - `work_clip`
   - `channel_step`
   - `work_step`

### 7.2 参数结构

涉及文件：

- `include/newcode/params.hpp`

需要做的事：

1. 删除或保留旧字段但标记兼容用途；
2. 新增 `CHANNEL_LLR_*` 与 `WORK_LLR_*` 两套字段；
3. 补公共步长推导辅助函数。

### 7.3 pipeline 层

涉及文件：

- `src/common/pipeline/pipeline_runner.cpp`
- `include/newcode/decoder_api.hpp`

需要做的事：

1. 只对 `channel_llr` 做动态 clip 计算；
2. 先得到 `channel_clip`；
3. 若开启“同步长 work”，则自动计算：

```text
work_clip = channel_clip * Q(work_bits) / Q(channel_bits)
```

4. 把两套 clip 和 bits 一起打包到 `DecodeRequest`。

### 7.4 decoder factory

第一阶段建议只支持：

- `chase_baseline`

涉及文件：

- `src/rx/decoder/factories/plain_decoder_factory.cpp`

建议新增一个新的解码入口，例如：

```cpp
decode_plain_split_qfloat<CHANNEL_BITS, WORK_BITS>(...)
```

它负责：

1. 用 `channel_bits/channel_clip` 量化 `request.channel_llr`
2. 初始化 `work_llr` 为 `WorkLLR`
3. 调用新的 mixed decoder 主流程
4. 输出 `post_decoder_llr` 时按 `work_clip` 反量化

### 7.5 OFEC 主流程模板

当前模板默认 `channel` 和 `work` 用同一个 `LLR` 类型。

涉及文件：

- `src/rx/ofec/detail/ofec_decode_impl.ipp`
- `src/rx/ofec/detail/ofec_window_impl.ipp`
- `src/rx/ofec/detail/ofec_tile_impl.ipp`
- `src/rx/ofec/detail/ofec_tile_input.ipp`
- `src/rx/ofec/detail/ofec_tile_writeback.ipp`

建议把模板拆成：

```cpp
template <typename ChannelLLR, typename WorkLLR>
```

其中：

- `channel_llr` 使用 `ChannelLLR`
- `work_llr` 使用 `WorkLLR`
- `tile_in` 使用 `WorkLLR`
- `ch_tile` 使用 `ChannelLLR`

### 7.6 `LinMatrixAdapter`

涉及文件：

- `include/newcode/ofec/common/lin_matrix_adapters.hpp`
- `src/rx/ofec/common/lin_matrix_adapters.ipp`

这是本次改造最关键的点。

当前接口默认：

```cpp
combine(const LLR& Lch, const LLR& La)
```

建议改成支持异类型：

```cpp
template <typename ChannelLLR, typename WorkLLR>
struct LinMatrixAdapter;
```

对“同步长分离量化”路径，推荐保留 code-domain 合成：

```cpp
return static_cast<float>(Lch.code() + La.code());
```

但要增加一个显式前提检查：

- 只有当 `channel_step == work_step` 时才允许这样做；
- 若以后支持“不同步长”的实验路径，则必须改为先变到真实 LLR 域再相加。

### 7.7 写回量化

涉及文件：

- `src/rx/ofec/detail/ofec_tile_decode.ipp`
- `src/rx/ofec/detail/ofec_tile_writeback.ipp`

需要做的事：

1. core 输出仍按 code-domain `float` 组织；
2. `ExtrinsicQuantizer` 要按 `WorkLLR` 的范围截断；
3. `tile_out` 写回的是 `WorkLLR`；
4. `last_tile_history_accum` 若继续保留 `float`，则要明确它保存的是“真实值”还是“code-domain 值”。

建议：

- `last_tile_history_accum` 继续保存真实 LLR 值；
- 这样最终输出和调试导出更直观。

## 8. 哪些地方仍然可以继续按 code 运算

只要满足“同步长”前提，下面这些地方都可以继续直接用 code：

1. `channel + work` 的主加法
2. `Lin` 的符号判定
3. `|Lin|` 的大小比较
4. least-reliable 位的排序
5. Chase 内部相关的度量输入

原因是：

- 同一步长时，`code` 与真实 LLR 只差一个公共比例因子；
- 对加法、比较、排序、判符号来说，这个公共比例因子不会改变结果。

## 9. 哪些地方需要重新检查

下面这些地方不应想当然复用旧逻辑：

1. 调试 CSV 中的“量化 code”导出
2. `post_decoder_llr` 的反量化口径
3. `work_llr` 导出的解释口径
4. `early-stop` 阈值是否是按真实值定义还是按 code-domain 定义
5. `ALPHA` 和 `beta` 的调参基准

特别注意：

- 当前有些逻辑名义上是 `float`，但在量化路径里其实承载的是 code-domain 值；
- 这次改造后，文档和变量命名最好把这一点说清楚，避免后续维护时混淆。

## 10. 建议的分阶段实施

### 阶段 1：只支持 `chase_baseline`

目标：

1. 跑通 `channel_bits != work_bits`
2. 保证 `channel_step == work_step`
3. 验证 BER 与 trace 行为

建议范围：

- 先只支持 `channel_bits = work_bits`
- 以及 `channel_bits = 6, work_bits = 7`

这样可以先把模板、参数和 factory 打通，再决定要不要做更一般的 pair 组合。

### 阶段 2：扩展到其他 decoder

后续再扩到：

- `topk_pruned`
- `global_pair`
- `group_minima`
- `ebchPF`

### 阶段 3：统一调试与导出

把下面这些导出统一成“channel 口径”和“work 口径”分开：

- channel quantized codes
- work quantized codes
- dequantized channel llr
- dequantized post-decoder llr

## 11. 验证建议

建议至少做以下检查。

### 11.1 参数检查

打印并确认：

```text
channel_bits = 6
work_bits = 7
channel_clip = 8
work_clip = 16.2580645
channel_step = work_step
```

### 11.2 单点数值检查

挑几个手工样本验证：

1. `channel_real = -8 -> channel_code = -31`
2. `work_real = +10 -> work_code ≈ +39`
3. `Lin_code = channel_code + work_code = +8`
4. `Lin_real ≈ 8 * delta`

### 11.3 回归检查

当：

- `channel_bits == work_bits`
- `work_clip == channel_clip`

时，新路径应尽量退化到旧行为。

### 11.4 BER 对比

建议至少比较三组：

1. 旧方案：`channel=6, work=6`
2. 新方案：`channel=6, work=7, same-step`
3. 浮点参考：`channel=16, work=16`

## 12. 本方案的结论

“保持 `work` 和 `channel` 的量化步长完全一样，只让 `work` 的可表示范围变大” 是一条合理且适合当前代码结构的改造方向。

它的核心优势是：

1. 保持量化精度不变；
2. 提升 `work` 的真实动态范围；
3. 允许继续沿用 code-domain 加法；
4. 改造范围可控；
5. 很适合先在 `chase_baseline` 上做第一版验证。

本方案最重要的实现原则只有一句：

```text
先固定公共步长 delta，再让 work 通过更大的 Q 扩真实范围。
```

对应你当前最关心的配置，就是：

- `channel: q6 @ clip = 8`
- `work: q7 @ clip = 8 * 63 / 31`

这会比简单设成 `work @ ±16` 更适合后续做 code-domain 加法与调试分析。
