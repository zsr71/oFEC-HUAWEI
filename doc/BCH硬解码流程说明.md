# BCH 硬解码流程说明

本文说明当前代码库里 BCH(255,239) 硬输入/硬输出译码的具体逻辑。

主要实现位置：

- `include/newcode/common/bch/bch_255_239.hpp`
- `src/common/bch/bch_255_239_decode.cpp`

主要接口：

```cpp
bool bch_255_239_decode_hiho_cw_255(
    const uint8_t* in255,
    uint8_t* out255,
    int* corrected_errors = nullptr,
    Bch255239DecodeTrace* trace = nullptr);
```

这里的 `hiho` 表示 hard-input / hard-output：

- 输入是 255 个硬判 bit，取值为 `0/1`
- 输出也是 255 个硬判 bit，表示 BCH 校正后的码字
- 这个接口只处理 BCH(255,239) 的前 255 位
- 扩展 BCH 的第 256 位 overall parity 不在这个函数内部处理

如果上层有 256 位码字，当前代码一般做法是：

- 前 255 位送入 BCH hard-decode
- 第 256 位 overall parity 由上层单独检查或重新计算

## 1. BCH 参数和数据约定

当前硬解码实现对应：

- BCH 码长：`255`
- 信息位：`239`
- BCH 校验位：`16`
- 纠错能力：`t = 2`
- 有限域：`GF(2^8)`
- primitive polynomial：`0x11D`

在 `Params` 里，扩展码长是 `256`：

- `BCH_N = 256`
- `BCH_K = 239`
- `BCH_OVERALL_IDX = 255`

但是 `bch_255_239_decode_hiho_cw_255(...)` 只处理前 255 位，不处理第 256 位。

上层把 LLR 转硬判时使用的约定通常是：

```cpp
LLR >= 0 -> bit 0
LLR <  0 -> bit 1
```

也就是说，正 LLR 表示更倾向于 0，负 LLR 表示更倾向于 1。

## 2. 有限域 GF(256) 初始化

实现内部有一个 `GF256` 结构，负责有限域运算。

使用的 primitive polynomial 是：

```cpp
0x11D
```

也就是：

```text
x^8 + x^4 + x^3 + x^2 + 1
```

初始化时会构造两个表：

- `alpha_to[i]`
  - 表示 `alpha^i` 对应的 GF(256) 元素
- `index_of[x]`
  - 表示 GF(256) 元素 `x` 对应的指数 `i`
  - 即 `x = alpha^i`

`0` 没有对数，所以：

```cpp
index_of[0] = -1;
```

有限域里的加法是异或：

```cpp
x + y = x ^ y
```

乘法通过指数表完成：

```cpp
x * y = alpha_to[(index_of[x] + index_of[y]) % 255]
```

除法通过指数相减完成：

```cpp
x / y = alpha_to[(index_of[x] - index_of[y] + 255) % 255]
```

## 3. Syndrome 计算

函数：

```cpp
compute_syndromes_1_4(const uint8_t* in255, uint8_t S[4])
```

会计算 4 个 syndrome：

```text
S1, S2, S3, S4
```

代码中对应：

```cpp
S[0] = S1
S[1] = S2
S[2] = S3
S[3] = S4
```

计算逻辑是遍历 255 个 bit。

如果第 `j` 位是 0，则跳过。

如果第 `j` 位是 1，则把对应的有限域元素累加到 syndrome 中：

```text
S1 ^= alpha^(1*j)
S2 ^= alpha^(2*j)
S3 ^= alpha^(3*j)
S4 ^= alpha^(4*j)
```

代码里写成：

```cpp
int e1 =  j       % 255;
int e2 = (j + e1) % 255;   // 2*j
int e3 = (j + e2) % 255;   // 3*j
int e4 = (j + e3) % 255;   // 4*j
```

然后：

```cpp
S[0] ^= alpha_to[e1];
S[1] ^= alpha_to[e2];
S[2] ^= alpha_to[e3];
S[3] ^= alpha_to[e4];
```

如果：

```text
S1 = S2 = S3 = S4 = 0
```

则说明当前 255 位已经满足 BCH 校验。

此时 hard-decode 会直接返回成功，不做纠错：

- `out255 = in255`
- `corrected_errors = 0`
- 返回 `true`

## 4. 主 hard-decode 流程

主函数流程如下：

```text
输入 in255
  |
  v
复制到本地 cw[255]
  |
  v
计算输入 syndrome S1..S4
  |
  +-- syndrome 全 0 -> 直接成功返回
  |
  v
Berlekamp-Massey 求错误定位多项式 sigma(x)
  |
  +-- L > 2 或异常 -> 失败返回
  |
  v
Chien search 找错误位置并翻转
  |
  +-- 找到的根数量不等于 L -> 失败返回
  |
  v
重新计算 syndrome S2
  |
  +-- 全 0 -> 成功返回
  |
  +-- 非 0 -> 失败返回
```

其中 `L` 是 Berlekamp-Massey 算出来的错误个数估计，也就是错误定位多项式的阶数。

由于当前 BCH 的纠错能力是 `t = 2`，所以只允许：

```text
L = 0, 1, 2
```

如果超过 2，当前实现不会尝试修正。

## 5. Berlekamp-Massey 部分

函数：

```cpp
berlekamp_massey_t2(const uint8_t S[4], uint8_t sigma[3])
```

作用是根据 `S1..S4` 求错误定位多项式：

```text
sigma(x) = 1 + c1*x + c2*x^2
```

代码里：

```cpp
sigma[0] = 1
sigma[1] = c1
sigma[2] = c2
```

核心变量：

- `C`
  - 当前错误定位多项式
- `B`
  - 上一次有效的辅助多项式
- `L`
  - 当前估计的错误个数
- `m`
  - 距离上一次更新的步数
- `b`
  - 上一次非零 discrepancy
- `d`
  - 当前 discrepancy

循环处理 `S1..S4`。

每一步先计算 discrepancy：

```text
d = S[n] + C1*S[n-1] + C2*S[n-2]
```

在 GF(256) 中，加法是异或，乘法是有限域乘法。

如果 `d == 0`：

- 当前多项式仍然匹配
- 不更新 `C`
- `m++`

如果 `d != 0`：

- 用 `d / b` 计算更新比例
- 根据 `m` 更新 `C`
- 如果 `2*L <= n`，说明需要提升错误定位多项式阶数
- 更新 `L`
- 保存旧的 `C` 到 `B`
- 更新 `b = d`
- 重置 `m = 1`

最后返回 `L`。

当前实现是 t=2 的轻量版本，只维护二阶多项式：

```cpp
uint8_t C[3] = {1,0,0};
uint8_t B[3] = {1,0,0};
```

所以它不是通用 BCH 任意 t 的 BM 实现，而是专门为 t=2 写的。

## 6. Chien Search 和纠错

函数：

```cpp
chien_and_correct(uint8_t* cw255, const uint8_t sigma[3], int L)
```

作用是：

1. 根据 `sigma(x)` 找错误位置
2. 找到后直接翻转 `cw255` 对应 bit

如果 `L == 0`：

- 表示 BM 认为没有错误
- 直接返回 `0`

否则从 `i = 1` 到 `255` 做 Chien search。

每一步计算：

```text
q = 1 + sigma1 * alpha^i + sigma2 * alpha^(2i)
```

如果：

```text
q == 0
```

说明这里找到一个错误定位多项式的根。

代码把根映射成 bit 位置：

```cpp
pos = 255 - i
```

注释里说明这个映射是为了和 AFF3CT 的位置定义一致。

找到所有错误位置后，要求：

```text
找到的根数量 == L
```

如果数量不一致：

- 返回 `-1`
- 表示 Chien search 失败

如果数量一致：

- 翻转对应位置：

```cpp
cw255[pos] ^= 1u;
```

- 返回实际修正数量

## 7. 纠错后的二次校验

Chien search 翻转 bit 后，主函数会再次计算 syndrome：

```cpp
compute_syndromes_1_4(cw, S2);
```

如果：

```text
S2[0] | S2[1] | S2[2] | S2[3] == 0
```

说明纠错后的码字通过 BCH 校验。

此时：

- `out255 = cw`
- `corrected_errors = corr`
- 返回 `true`

其中 `corr` 是 Chien search 实际翻转的 bit 数。

如果二次 syndrome 仍然非 0：

- `out255 = cw`
- `corrected_errors = -1`
- 返回 `false`

也就是说，即使返回失败，`out255` 也可能包含一次尝试纠错后的结果。

上层一般应该以返回值 `true/false` 为准，而不是只看 `out255`。

## 8. 失败分支和返回值

主接口的返回语义如下。

| 场景 | 返回值 | `corrected_errors` | `out255` |
|---|---:|---:|---|
| 输入 syndrome 全 0 | `true` | `0` | 原始输入 |
| BM 求出的 `L` 不合法 | `false` | `-1` | 原始输入副本 |
| Chien search 根数量不匹配 | `false` | `-1` | 原始输入 |
| 纠错后二次 syndrome 非 0 | `false` | `-1` | 尝试纠错后的 `cw` |
| 成功纠错 | `true` | `1` 或 `2` | 修正后的码字 |

需要注意：

- 当前纠错能力最多 2 bit
- 超过 2 bit 错误时，可能返回失败
- 某些超过 2 bit 的错误模式也可能被误校正，这属于 BCH 硬译码本身的风险，需要靠外层 metric、overall parity 或 Chase 竞争来降低影响

## 9. Trace 输出

现在接口支持一个可选 trace：

```cpp
struct Bch255239DecodeTrace {
    std::array<uint8_t,4> input_syndromes{};
    std::array<uint8_t,4> output_syndromes{};
    bool has_output_syndromes = false;
};
```

如果调用时传入：

```cpp
trace != nullptr
```

则 hard-decode 会记录：

- `input_syndromes`
  - 输入码字进入 BCH 前的 `S1..S4`
- `output_syndromes`
  - BCH 尝试纠错后的 `S1..S4`
- `has_output_syndromes`
  - 表示是否已经写入输出 syndrome

如果输入本来 syndrome 全 0：

- `output_syndromes = input_syndromes`
- `has_output_syndromes = true`

如果 BM 或 Chien search 早期失败：

- 可能只有 `input_syndromes`
- `has_output_syndromes = false`

这个 trace 当前主要用于 Chase candidate syndrome CSV：

- 每个 Chase 测试样本都会先过 BCH hard-decode
- BCH 内部本来就计算了 `S1..S4`
- CSV 只取 `input_syndromes[0]` 和 `input_syndromes[2]`
- 也就是保存 `S1` 和 `S3`

正常不导出 CSV 时，调用方传 `nullptr`，不会保存这些 trace。

## 10. 辅助检查函数

除了 hard-decode，当前还提供了几个 syndrome 检查接口。

### 10.1 `bch_255_239_syndromes_zero_cw_255`

```cpp
bool bch_255_239_syndromes_zero_cw_255(const uint8_t* in255)
```

作用：

- 只计算 `S1..S4`
- 判断是否全 0
- 不做纠错

用途：

- early-stop 条件判断
- TPC 行列校验

### 10.2 `bch_255_239_syndromes_1_4_cw_255`

```cpp
std::array<uint8_t,4> bch_255_239_syndromes_1_4_cw_255(const uint8_t* in255)
```

作用：

- 返回完整 `S1..S4`

用途：

- 调试
- 统计
- 需要 syndrome 原始 GF(256) 数值的分析

### 10.3 `bch_255_239_syndrome_nonzero_mask_cw_255`

```cpp
uint8_t bch_255_239_syndrome_nonzero_mask_cw_255(const uint8_t* in255)
```

作用：

- 把 `S1..S4` 是否非 0 压成一个 4-bit mask

bit 含义：

```text
bit0 = S1 != 0
bit1 = S2 != 0
bit2 = S3 != 0
bit3 = S4 != 0
```

这个 mask 不包含 syndrome 的具体数值，只表示每个 syndrome 分量是否为 0。

它适合用于：

- early-stop 细节统计
- MUX 优先级调度里的“校验式坏得多不多”这类粗粒度判断

它不适合用于：

- `S1^3 == S3` 这种需要 syndrome 原始数值的分析

## 11. 当前主要调用场景

### 11.1 Chase 软解码候选

位置示例：

- `src/rx/ofec/plain/detail/chase256_plain_impl.ipp`
- `src/rx/ofec/plain/detail/chase256_topk_pruned_impl.ipp`
- `src/rx/ofec/plain/detail/chase256_global_pair_impl.ipp`
- `src/rx/ofec/plain/detail/chase256_group_minima_impl.ipp`
- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp`

Chase 会生成 `CHASE_NTEST` 个测试样本。

每个测试样本流程是：

```text
当前硬判 hard_ch
  |
  v
按 test pattern 翻转最不可靠位
  |
  v
前 255 位送入 BCH hard-decode
  |
  v
得到候选 BCH 码字 cw255
  |
  v
补 overall parity 得到 256 位候选
  |
  v
计算 metric，参与 ML / competing codeword 选择
```

BCH hard-decode 在这里的作用是把 Chase 生成的测试样本拉回到合法 BCH 码字附近。

### 11.2 early-stop action mode3

位置：

- `src/rx/ofec/earlystop/row_early_stop_process_3.ipp`

mode3 流程是：

```text
lin256 硬判成 256 bit
  |
  v
前 255 位 BCH hard-decode
  |
  +-- BCH 失败 -> mode3 返回 false
  |
  v
重新计算 overall parity
  |
  v
输出固定幅度 ±EARLY_STOP_ACTION_HARD_LLR_MAG
```

这里 BCH hard-decode 失败会导致 mode3 动作失败。

这和 mode6 不同：

- mode6 不做 BCH hard-decode
- mode6 直接按 `sign(lin)` 输出 `±EARLY_STOP_ACTION_SIGN_BETA`

### 11.3 OFEC hard decode fallback

位置：

- `src/rx/ofec/ofec_hard_decode_bch.cpp`

流程是：

```text
Lin256 硬判成 256 bit
  |
  v
前 255 位 BCH hard-decode
  |
  +-- 失败 -> 返回 false
  |
  v
补第 256 位 overall parity
  |
  v
用固定幅度 HARD_LLR_MAG 构造后验 Lpost
  |
  v
Y2 = Lpost - Lin256
```

这里输出的 `Y2` 是一种硬判回退外信息。

### 11.4 TPC 校验

位置：

- `src/tpc/tpc_decoder.cpp`

TPC 里使用的是 syndrome 检查，不是纠错：

```text
前 255 位 BCH syndrome 全 0
并且
第 256 位 overall parity 正确
```

二者都满足，才认为这一行或这一列通过校验。

## 12. 当前实现需要注意的点

### 12.1 BCH hard-decode 不管 overall parity

`bch_255_239_decode_hiho_cw_255(...)` 只处理 255 位 BCH。

第 256 位 overall parity 必须由上层处理。

因此，如果要判断完整 256 位扩展 BCH 是否通过，不能只看 BCH syndrome，还要检查 overall parity。

### 12.2 返回失败时不要信任 `out255`

失败场景下，`out255` 可能是：

- 原始输入
- 原始输入副本
- 尝试纠错后的中间结果

所以调用方应该先看返回值。

只有返回 `true` 时，`out255` 才能作为成功译码结果使用。

### 12.3 当前是 t=2 专用实现

BM 和 Chien search 都是围绕 `t=2` 写的。

如果以后换 BCH 参数或纠错能力，这里不能直接复用。

需要重新检查：

- syndrome 数量
- BM 多项式阶数
- Chien search 逻辑
- corrected error 上限

### 12.4 syndrome 数值是 GF(256) 元素

`S1..S4` 的类型是 `uint8_t`，但它们不是普通整数意义下的 syndrome。

它们表示 GF(256) 元素。

所以在 Matlab 里分析 `S1^3 == S3` 时，需要按 GF(256) 乘法计算，不应该用普通整数三次方。

