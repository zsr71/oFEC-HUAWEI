# MUX优先级规则说明

## 背景
`kMuxPriorityRule` 只在下面这个条件同时成立时才生效：

- `kMuxSchedulingMode = 1`

这时，MUX 不再按原始顺序裁掉 `NeedSiso` 的码字，而是会先根据 early-stop 细节给每个 `NeedSiso row` 打分，再按分数排序。

如果：

- `kMuxSchedulingMode = 0`

那么 `kMuxPriorityRule` 不生效，系统仍然沿用旧的顺序裁剪逻辑。

## 参与排序的对象
当前排序只作用于：

- `NeedSiso`

不会作用于：

- `EarlyStopped`
- `Unscheduled`

也就是说：

- 已经命中 early-stop 的 row，不参与这轮优先级竞争
- 真正参加排序的是“还需要正常 SISO/Chase”的那些 row

## 排序输入信息
当前优先级打分使用的是每个 row 在 early-stop 条件 1 检测后留下的细节信息：

- `syndrome_bits`
- `bch_passed`
- `overall_parity_passed`

其中：

- `syndrome_bits`
  - 是 BCH syndrome 的原始 bit 信息
  - 当前代码里实际使用的是其中有多少个非 0

- `bch_passed`
  - 表示 BCH syndrome 是否全 0

- `overall_parity_passed`
  - 表示 overall parity 是否通过

## 当前实现位置
优先级打分与排序实现在：

- [mux_siso_budget.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_siso_budget.cpp)
- [mux_group_budget.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_group_budget.cpp)

其中核心函数是：

- `early_stop_priority_score(...)`
- `trim_need_indices_by_priority(...)`

## 规则 0：`kMuxPriorityRule = 0`
这个规则的含义是：

- 更差的 row 优先

当前分数公式是：

```cpp
score = syndrome_nonzero_count * 10 + (overall_parity_passed ? 0 : 1)
```

其中：

- `syndrome_nonzero_count`
  - 表示 syndrome 中非 0 的数量

- `overall_parity_passed ? 0 : 1`
  - 如果 overall parity 没过，再额外加 1 分

### 直观解释
这个规则认为：

- syndrome 非 0 越多，说明这行离 early-stop 条件 1 越远
- overall parity 也不过，则再稍微更“差”一点

所以最终会优先保留：

- syndrome 更差的 row

### 例子
假设有两行：

- Row A
  - `syndrome_nonzero_count = 5`
  - `overall_parity_passed = false`
  - `score = 5 * 10 + 1 = 51`

- Row B
  - `syndrome_nonzero_count = 2`
  - `overall_parity_passed = true`
  - `score = 2 * 10 + 0 = 20`

那么：

- Row A 优先级更高

## 规则 1：`kMuxPriorityRule = 1`
这个规则的含义是：

- 更接近通过 early-stop 的 row 优先

当前分数公式是：

```cpp
if (bch_passed && !overall_parity_passed) {
  score = 100;
} else {
  score = (bch_passed ? 50 : 0) +
          (overall_parity_passed ? 5 : 0) -
          syndrome_nonzero_count;
}
```

### 直观解释
这个规则的判断顺序是：

1. 如果 BCH 已经通过，只差 overall parity
   - 直接给最高优先级 `100`

2. 否则：
   - BCH 通过会加很多分
   - overall parity 通过会再加一点分
   - syndrome 非 0 越多，反而扣分

所以最终会优先保留：

- 已经很接近通过 early-stop 的 row

### 例子
假设有两行：

- Row A
  - `bch_passed = true`
  - `overall_parity_passed = false`
  - `score = 100`

- Row B
  - `bch_passed = false`
  - `overall_parity_passed = true`
  - `syndrome_nonzero_count = 2`
  - `score = 0 + 5 - 2 = 3`

那么：

- Row A 优先级明显更高

## 排序方式
排序使用的是：

- `std::stable_sort`

并且按：

- 分数从高到低

排序。

这意味着：

- 分数高的 row 会排在前面
- 如果两个 row 分数相同，则保持原始顺序不变

## 排序后怎么裁剪
排序只决定“谁排前面”，真正的预算裁剪仍然是：

- 保留前 `keep_budget` 个
- 后面的改成 `Unscheduled`

例如：

- 某个 tile 里一共有 12 个 `NeedSiso`
- 当前只允许保留 8 个

那么流程就是：

1. 先给这 12 个 row 打分
2. 按分数从高到低排序
3. 保留前 8 个
4. 后 4 个改成 `Unscheduled`

## `group_g > 1` 时的逻辑
如果：

- `kMuxGroupG > 1`

那么不是全 tile 一起排序，而是：

1. 先把全部 row 按 group 分组
2. 再把预算均分到各组
3. 每个组内部各自按 `kMuxPriorityRule` 排序
4. 每组内部再裁掉超预算的 `NeedSiso`

也就是说：

- `kMuxPriorityRule` 仍然生效
- 但生效范围变成“组内排序”

## 和 early-stop 的关系
`kMuxPriorityRule` 并不会改变：

- 哪些 row 已经 `EarlyStopped`

它只影响：

- 那些还没 early-stop、仍属于 `NeedSiso` 的 row，在预算不够时谁优先保留

所以它的作用可以概括成：

- early-stop 决定谁直接放行
- `kMuxPriorityRule` 决定剩下需要 SISO 的那些 row 谁优先拿到解码资源

## 建议理解
可以直接这样记：

- `kMuxPriorityRule = 0`
  - 更差的先解

- `kMuxPriorityRule = 1`
  - 更接近通过 early-stop 的先解

## 一句话总结
`kMuxPriorityRule` 是 `kMuxSchedulingMode = 1` 时使用的二级规则，用来决定：

- 在预算不足时，哪些 `NeedSiso` row 应该优先被保留进入正常 SISO/Chase 解码。  
