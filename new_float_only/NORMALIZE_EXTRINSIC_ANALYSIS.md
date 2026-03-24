# `normalize_extrinsic` 与 MATLAB TPC 归一化对比分析

本文分析 [`new_float_only`](/home/zsr71/projects/newcode/new_float_only) 里的 `normalize_extrinsic` 开关，和用户给出的 MATLAB `decode()` / `compBlkDecode()` 实现中的归一化是否一致。

结论先写在前面：

- 两者的**目标相似**：都是为了把外信息 `extrinsic` 的平均幅度控制在稳定范围，避免迭代过程中幅度发散或过弱。
- 两者的**实现层级不一样**：MATLAB 是在每次行/列半轮结束后，对整块 `w` 做归一化；`new_float_only` 是在每个 tile 的 Chase 输出 `lout` 上做局部归一化。
- 两者**不是完全相同的实现**，只能说是“同类思路，但不是同一公式的逐行复现”。

## 1. MATLAB 代码里的归一化做了什么

用户给出的 MATLAB 代码里，归一化其实有两层。

### 1.1 输入 LLR 归一化

在 `decode()` 一开始：

```matlab
llrMean = mean(abs(llrMat(:)));
if(llrMean ~= 0)
    llrMat = llrMat/llrMean;   % LLR normalization
end
```

这一步的含义是：

- 对整个输入信道 LLR 矩阵 `llrMat` 取平均绝对值
- 如果均值非零，就把整个 `llrMat` 除以这个均值

所以 MATLAB 版本会先把**信道输入**的平均幅度归一化到大约 1。

### 1.2 每个半轮后的外信息归一化

在行译码结束后：

```matlab
mean_abs_w = mean(abs(w(abs(w) ~= beta(decstpidx))));
if(~isnan(mean_abs_w) && mean_abs_w ~= 0)
    w = w/mean_abs_w;
end
```

列译码结束后也有完全相同的一段。

这一步的含义是：

1. 当前半轮得到整块外信息矩阵 `w`
2. 把那些数值**恰好等于** `beta(decstpidx)` 的元素排除掉
3. 对剩余元素的绝对值求平均
4. 用这个均值去归一化整个 `w`

所以 MATLAB 的做法可以概括成：

- 先对输入 `llrMat` 做全局归一化
- 再对每次半轮后的整个 `w` 做全局归一化

## 2. `new_float_only` 里的归一化做了什么

`new_float_only` 的开关来自：

- [`pipeline_runner.cpp`](/home/zsr71/projects/newcode/new_float_only/src/pipeline/pipeline_runner.cpp)
- [`ofec_frame_decode.cpp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/ofec_frame_decode.cpp)
- [`ofec_tile_decode.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_decode.ipp)

调用链是：

1. app 层设置 `config.normalize_extrinsic`
2. `run_pipeline()` 把它传给 `decode_plain_llr(...)`
3. 每个 tile 的 `decode_tile(...)` 决定是否做外信息归一化

### 2.1 归一化只作用于 tile 的 `lout`

在 [`ofec_tile_decode.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_decode.ipp) 中：

```cpp
if (normalize_extrinsic && !use_hard_decode)
{
  normalize_extrinsic_lout(decoder_res.lout, decoder_res.produced_rows, p.beta);
}
```

可以看到：

- 归一化对象不是整块矩阵 `w`
- 而是当前 tile 内部当前一次 Chase decoder 输出的 `lout`
- 并且只有 soft decode 分支会做，hard decode 不做

### 2.2 归一化公式

实际函数在 [`ofec_tile_extrinsic_normalize.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_extrinsic_normalize.ipp)：

```cpp
double acc = 0.0;
std::size_t cnt = 0;
...
acc += ...
++cnt;
...
const Float g_alpha = static_cast<Float>(acc / static_cast<double>(cnt));
const Float scale = static_cast<Float>(1.0f) / g_alpha;
lout[r][j] *= scale;
```

所以核心也是：

```text
scale = 1 / mean_abs
lout = lout * scale
```

也就是把当前这批外信息的平均绝对值拉回到大约 1。

### 2.3 对 `beta` fallback 的处理方式

这里和 MATLAB 最重要的不同点在于：`new_float_only` 没有简单地“排除所有等于 beta 的值”，而是做了一个近似判断：

```cpp
is_fallback(w) = fabs(fabs(w) - beta) <= tol
```

如果某个值被认为是 fallback 值，那么统计均值时不是加 `abs(w)`，而是加：

```cpp
abs(w / beta)
```

如果 `w = ±beta`，那它的贡献就是 `1`。

所以这一步的效果更像：

- 普通外信息：按真实幅度计入均值
- fallback 外信息：按“单位 fallback”计入均值

这和 MATLAB 的“直接排除 `abs(w) == beta` 的元素”并不相同。

## 3. 归一化之后还做了什么

在 [`ofec_tile_decode.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_tile_decode.ipp) 中，归一化后还有一步：

```cpp
decoder_res.lout[r][j] *= p.ALPHA;
```

也就是说，`new_float_only` 当前的顺序是：

1. 先对 tile 的 `lout` 做平均幅度归一化
2. 再乘当前 tile 的 `ALPHA`

而且这个 `ALPHA` / `beta` 是每个 tile 单独取的，来自：

- [`ofec_window_impl.ipp`](/home/zsr71/projects/newcode/new_float_only/src/rx/ofec/detail/ofec_window_impl.ipp)

```cpp
tile_params.beta = pick_float(p.beta_list, t, p.beta);
tile_params.ALPHA = pick_float(p.ALPHA_LIST, t, p.ALPHA);
```

所以当前实现的实际形式更接近：

```text
extrinsic_out(tile t) = ALPHA_t * normalize(extrinsic_raw(tile t), beta_t)
```

## 4. 两者相同点

两者相同点主要有这些：

1. 都是在控制外信息幅度
2. 都是基于“平均绝对值”做缩放
3. 都和 Pyndiah/Chase 风格的 `beta` fallback 机制有关
4. 归一化后的目标都可以理解为“把典型外信息幅度拉回到 1 附近”

所以如果只看算法意图，它们是相近的。

## 5. 两者不同点

两者不同点更关键：

### 5.1 归一化层级不同

MATLAB：

- 先归一化整块输入 `llrMat`
- 再在每个行/列半轮后归一化整块 `w`

`new_float_only`：

- 不对整个输入 `channel_llr` 做那种 `mean(abs(llrMat(:)))` 归一化
- 只在每个 tile 的 `lout` 上做局部归一化

### 5.2 `beta` 元素处理不同

MATLAB：

- 直接排除 `abs(w) == beta` 的元素

`new_float_only`：

- 用容差判断“是否接近 beta”
- 若接近，则按 `abs(w / beta)` 计入均值

所以这部分不是同一实现。

### 5.3 `alpha` 的使用方式不同

MATLAB 代码里，`alpha(decstpidx)` 是在 `softin = r + alpha*win` 中参与下一步软输入合成。

`new_float_only` 当前的 `ALPHA` 则是在 tile 外信息归一化之后，再直接对 `lout` 做缩放：

```cpp
lout *= p.ALPHA;
```

所以虽然名字相同，但放置的位置并不完全一致。

### 5.4 调度结构不同

MATLAB 版本是标准 TPC：

- 整块行译码
- 整块归一化
- 整块列译码
- 整块归一化

`new_float_only` 当前实现是：

- 滑窗
- tile
- tile 内部 row-level Chase core
- tile 级外信息写回

因此即使公式接近，归一化发生的时空位置也不同。

## 6. 能不能说“和 MATLAB 一样”

不能严格说“一样”。

更准确的说法应该是：

- `new_float_only` 的 `normalize_extrinsic` 和 MATLAB TPC 代码**属于同类归一化思想**
- 目标都是把 Pyndiah 外信息的平均幅度控制在 1 附近
- 但它**不是 MATLAB 实现的逐步复刻**

如果要用一句最短的话概括：

> 思路相近，层级不同，公式也不完全相同。

## 7. 如果想更接近 MATLAB 版本，需要改哪些点

如果后面希望 `new_float_only` 更贴近 MATLAB 这段代码，主要要改 3 类地方：

1. 在进入解码前，对整块 `channel_llr` 先做一次 `mean(abs(.))` 全局归一化
2. 不要在 tile 级做归一化，而改成在一个更高层的“完整半轮输出”上归一化
3. `beta` 处理改成更接近 MATLAB 的：
   - 直接排除 `abs(w) == beta` 的元素
   - 再对剩余元素求均值

但要注意，一旦这样改，行为就会明显更接近“整块二维 TPC 迭代器”，而不是当前这个滑窗/tile 版本。

## 8. 当前建议理解

对于当前 `new_float_only`，建议把 `kNormalizeExtrinsic` 理解成：

- 一个**tile 级的外信息平均幅度归一化开关**
- 用于让不同 tile、不同 pattern 下的外信息尺度更稳定
- 它参考了 Pyndiah 风格的思路，但不是 MATLAB 那份 TPC 代码的逐句等价实现
