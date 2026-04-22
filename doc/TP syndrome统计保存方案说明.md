# TP / Chase 候选 syndrome 统计保存方案说明

本文说明一个问题：

- 对于**进入 Chase 软解码的非 early-stop 码字**
- 每个码字会生成 `CHASE_NTEST` 个测试样本
- 现在希望把这些**候选测试样本**的 `S1`、`S3` 记录下来
- 后面在 Matlab 里统计：
  - `s1 == 0`
  - `s1^3 == s3`
  - 两者同时满足

这里讨论的对象是 **Chase 候选码字级别**，不是 early-stop row 级别。

## 1. 需要记录的对象

当前 Chase 软解码会对每个非早停 row：

- 枚举 `CHASE_NTEST` 个测试样本
- 对每个测试样本得到一个候选码字

你要记录的是：

- **每个候选测试样本的 syndrome**
- 至少保存：
  - `S1`
  - `S3`

这样 Matlab 里就可以后处理：

- `s1_zero = (s1 == 0)`
- `s1_cubed_eq_s3 = gf_pow(s1, 3) == s3`
- `both = s1_zero && s1_cubed_eq_s3`

## 2. 为什么不是只存 `syndrome_bits`

当前程序里的 `syndrome_bits` 是 early-stop 层的一个**非零掩码**：

- bit0: `S1 != 0`
- bit1: `S2 != 0`
- bit2: `S3 != 0`
- bit3: `S4 != 0`

它记录的是：

- 哪几个 syndrome 分量为 0
- 哪几个 syndrome 分量非 0

它**不记录** syndrome 的原始 GF(256) 数值。

所以：

- `syndrome_bits` 对 early-stop row 的粗粒度统计有用
- 但对 `s1^3 == s3` 这种数值关系**不够**

## 3. 这个统计能做到什么粒度

如果额外保存 `S1` 和 `S3`，那么后面可以做到这些粒度：

### 3.1 候选样本粒度

对每一个 Chase 测试样本单独记录：

- `candidate_idx`
- `S1`
- `S3`
- `s1_zero`
- `s1_cubed_eq_s3`

### 3.2 码字粒度

对一个非 early-stop codeword 的所有 `CHASE_NTEST` 候选做统计：

- `S1 == 0` 的候选数
- `S1^3 == S3` 的候选数
- 两者同时满足的候选数

### 3.3 tile / invocation 粒度

如果 CSV 里同时带上：

- `ebn0_db`
- `seed_index`
- `invocation`
- `tile_index`
- `row / bit target`

那 Matlab 里还可以进一步按：

- tile
- invocation
- seed

做聚合统计。

## 4. 当前代码里应该在哪里记录

这个统计应该放在 **Chase 候选生成阶段**，不是 early-stop 判定阶段。

最相关的代码位置是：

- `src/rx/ofec/plain/detail/chase256_plain_impl.ipp`
- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp`
- 以及其它 Chase 变体对应的实现文件

原因是：

- 这里本来就有 `for (int c = 0; c < NTEST; ++c)`
- 每个 `c` 就是一个测试样本
- 每个测试样本都会调用 BCH hard-decode
- BCH hard-decode 内部本来就已经计算了 syndrome

当前更准确的做法不是在 Chase 外面额外再算一遍 syndrome，而是：

1. 扩展 BCH hard-decode 接口
2. 把内部已经算出的 `S1..S4` 通过 trace 参数带出来
3. Chase 候选循环只负责保存这个 trace 结果

这样可以避免对每个候选重复计算 syndrome。

## 5. 推荐的接口方案

建议新增一个可选 trace 结构，例如：

```cpp
struct Bch255239DecodeTrace {
  std::array<uint8_t, 4> input_syndromes{};
  std::array<uint8_t, 4> output_syndromes{};
  bool has_output_syndromes = false;
};
```

然后给 BCH hard-decode 增加一个可选参数：

```cpp
bool bch_255_239_decode_hiho_cw_255(
    const uint8_t* in255,
    uint8_t* out255,
    int* corrected_errors,
    Bch255239DecodeTrace* trace);
```

其中：

- `input_syndromes`
  - 对应测试样本进入 BCH 前的 syndrome
  - 这是当前 TP 统计最关心的对象
- `output_syndromes`
  - 对应 BCH 修正后的候选码字 syndrome
  - 主要用于调试 BCH hard-decode 是否成功
- `trace == nullptr`
  - 表示不记录 syndrome
  - 正常解码路径继续走原逻辑

第一版建议 Chase 侧记录：

- `input_syndromes[0]`，也就是 `S1`
- `input_syndromes[2]`，也就是 `S3`

## 6. 性能影响控制

这个方案的关键目标是：**记录开启时拿到 `S1/S3`，记录关闭时不影响解码性能。**

为此建议这样设计：

- `trace` 参数默认为 `nullptr`
- 只有打开 candidate syndrome dump 时，Chase 才传入有效 trace 指针
- BCH decode 内部本来就会计算 `input_syndromes`
- 因此打开记录时只是多做一次内存拷贝，不额外计算 syndrome
- 关闭记录时只多一个空指针判断

也就是说：

- 正常 BER sweep / single 跑性能时，保持 `trace == nullptr`
- 只有调试导出 CSV 时，才记录 `S1/S3`

这样不会改变主解码结果，也不会对普通运行引入明显性能开销。

## 7. 推荐的保存方式

第一版最小增量建议：

- 保留当前的 `syndrome_bits`
- 额外保存 `S1`
- 额外保存 `S3`

更完整一点也可以直接保存：

- `S1`
- `S2`
- `S3`
- `S4`

但如果只是为了后面 Matlab 分析 `s1 == 0` 和 `s1^3 == s3`，第一版保存 `S1 + S3` 就够了。

## 8. 方案对比

### 方案 A：只存 `syndrome_bits`

优点：

- 体积小
- 改动少

缺点：

- 只能做零 / 非零模式统计
- 不能精确算 `s1^3 == s3`

适合：

- 只看 syndrome 是否通过
- 只看粗粒度异常分类

### 方案 B：存 `syndrome_bits + S1 + S3`

优点：

- 兼容当前粗粒度统计
- 又能做 `s1^3 == s3`

缺点：

- 比方案 A 多两个字段

适合：

- 既想保留当前实现
- 又想做候选级 syndrome 分析

### 方案 C：直接存完整 `S1..S4`

优点：

- 最灵活
- 后面想继续扩展分析，不用再改数据结构

缺点：

- CSV 更大
- 比方案 B 稍重

适合：

- 后续可能继续扩展 syndrome 相关实验

## 9. 如果要修改代码，涉及哪些文件

### 9.1 第一版建议修改的地方

1. `include/newcode/common/bch/bch_255_239.hpp`
   - 新增 `Bch255239DecodeTrace`
   - 新增带 trace 参数的 BCH decode 重载

2. `src/common/bch/bch_255_239_decode.cpp`
   - 在已有 `compute_syndromes_1_4(...)` 之后，把 `S1..S4` 写入 trace
   - 不改变原有 decode 判定逻辑

3. `src/rx/ofec/plain/detail/chase256_plain_impl.ipp`
   - 在 candidate loop 里按需传入 trace
   - 保存每个 candidate 的 `S1/S3`

4. 其它 plain Chase 变体
   - `src/rx/ofec/plain/detail/chase256_topk_pruned_impl.ipp`
   - `src/rx/ofec/plain/detail/chase256_global_pair_impl.ipp`
   - `src/rx/ofec/plain/detail/chase256_group_minima_impl.ipp`
   - 使用同一套 trace 接口保存每个 candidate 的 `S1/S3`

5. `src/rx/ofec/ebchPF/chase256_ebchPF.cpp`
   - 如果该 decoder 也需要同样统计，同步接入 trace

6. `src/rx/ofec/plain/detail/chase256_plain_trace_csv.ipp`
   - 把 candidate 级 `S1/S3` 导出到 CSV

### 9.2 如果要导出到 CSV

还需要改对应的 trace 导出逻辑，比如：

- `src/rx/ofec/plain/detail/chase256_plain_trace_csv.ipp`
- `src/rx/ofec/ebchPF/chase256_ebchPF.cpp`

给每个 candidate CSV 增加列：

- `S1`
- `S3`

第一版采用追加 candidate 明细表的方式，不新增大量控制台打印。

如果以后要做汇总分析，还可以再加：

- `s1_zero`
- `s1_cubed_eq_s3`

### 9.3 不建议第一版就动的地方

- `src/rx/ofec/earlystop/tile_early_stop_stats.ipp`
- `include/newcode/ofec/earlystop/tile_early_stop_result.hpp`
- `src/rx/ofec/mux/mux_siso_budget.cpp`
- `src/rx/ofec/mux/mux_group_budget.cpp`

原因是：

- 你现在要的是 **Chase 候选级统计**
- 不是 early-stop row 统计
- 也不是 MUX priority 逻辑

先把候选样本的 `S1/S3` 保存下来最稳。

## 10. 结论

如果你的目标是：

- 对每个进入 Chase 的非 early-stop 码字
- 记录 `CHASE_NTEST` 个测试样本的 `S1 / S3`
- 后面在 Matlab 里统计：
  - `s1 == 0`
  - `s1^3 == s3`
  - 两者同时满足

那当前最稳妥的方案是：

- **保留现有 `syndrome_bits`**
- **扩展 BCH decode trace，把已经算出的 `S1/S3` 带出来**
- **只在打开候选 syndrome dump 时记录**

这样既不破坏现有早停 / MUX 逻辑，也能满足后续分析需要。
