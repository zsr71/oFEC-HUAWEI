# `kChaseTopkKeep` 在解码流程中的生效逻辑

这份说明只针对 `decoder_name = "chase_topk_pruned"` 时，`kChaseTopkKeep` 是怎么从顶层参数一路传到实际解码逻辑里的。

## 1. 顶层参数入口

### `ofec_single`

在单次运行里，顶层常量是：

- [apps/ofec_single.cpp](/home/zsr71/projects/newcode/apps/ofec_single.cpp)
  - `kChaseTopkKeep`

它被写入：

- [include/newcode/ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp)
  - `ofec_single::Config::chase_topk_keep`

然后在参数构造阶段写进运行参数：

- [src/ofec_single/ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp)
  - `params.CHASE_TOPK_KEEP = cfg.chase_topk_keep;`

同时这里会做基本校验：

- `chase_topk_keep < 1` 会直接报错

对应字段定义在：

- [include/newcode/params.hpp](/home/zsr71/projects/newcode/include/newcode/params.hpp)
  - `int CHASE_TOPK_KEEP = 8;`

也就是说，`ofec_single` 下 `kChaseTopkKeep` 就是这次运行固定使用的 Top-K 保留数。

### `ofec_sweep`

在 sweep 里有两套入口：

- 固定值
  - [apps/ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)
    - `kChaseTopkKeep`
- 扫描候选
  - [apps/ofec_sweep.cpp](/home/zsr71/projects/newcode/apps/ofec_sweep.cpp)
    - `kChaseTopkKeepCandidates`

它们会进入：

- [src/ofec_sweep/ofec_sweep_scenarios.cpp](/home/zsr71/projects/newcode/src/ofec_sweep/ofec_sweep_scenarios.cpp)
  - `choose_candidates(...)`

这里的规则是：

- 如果 `kChaseTopkKeepCandidates` 非空，就优先用候选列表展开场景
- 如果 `kChaseTopkKeepCandidates` 为空，才回退使用单值 `kChaseTopkKeep`

所以 sweep 下的优先级是：

1. `kChaseTopkKeepCandidates`
2. `kChaseTopkKeep`

注意这组参数只在 `decoder_name == "chase_topk_pruned"` 时才真正展开；对其它 decoder，只会保留一个占位值，不参与场景扩展。

## 2. 什么时候这个参数会真正生效

`CHASE_TOPK_KEEP` 只有在：

- `decoder_name = "chase_topk_pruned"`

时才有实际作用。

这一点在 sweep 场景构造里也写死了：

- [src/ofec_sweep/ofec_sweep_scenarios.cpp](/home/zsr71/projects/newcode/src/ofec_sweep/ofec_sweep_scenarios.cpp)

逻辑是：

- 如果当前 decoder 是 `chase_topk_pruned`
  - `topk_keep_values = choose_candidates(config.chase_topk_keep_candidates, config.chase_topk_keep)`
- 否则
  - 只给一个固定值，不展开扫描

也就是说：

- `chase_baseline` 不看它
- `chase_global_pair` 不看它
- `chase_group_minima` 不看它
- `chase_overall_parity_search` 不看它

## 3. 进入 `chase_topk_pruned` 之后，Top-K 是怎么选的

真正的实现位于：

- [src/rx/ofec/plain/detail/chase256_topk_pruned_impl.ipp](/home/zsr71/projects/newcode/src/rx/ofec/plain/detail/chase256_topk_pruned_impl.ipp)

### 第一步：先生成全部候选

代码先按当前 `CHASE_L` 和 `CHASE_NTEST`：

- 找最不可靠位
- 生成全部 test pattern
- 得到 `NTEST` 个候选
- 每个候选都跑一次 BCH 硬解
- 计算该候选对应的 `score = -dist`

这里的 `comps` 里装的是全部候选，每个候选有：

- `score`
- `idx`
- `good`
- `corrected_errors`

其中：

- `good = true` 表示这个候选在 BCH 硬解码时成功
- `good = false` 表示这个候选虽然有分数，但没有通过 BCH 合法性筛选

### 第二步：先筛 `good`

当前代码不是在“全部候选”里直接取 Top-K，而是：

1. 先从 `comps` 里筛出 `good == true` 的候选
2. 得到 `kept_comps`

也就是说，现在 Top-K 的池子是：

- 仅 `good` 候选

不是：

- 全部候选

### 第三步：在 `good` 候选里按分数排序

排序规则是：

- 先按 `score` 降序
- 若分数相同，再按 `idx` 升序

对应含义：

- 分数更高的候选优先
- 分数相同就按原始候选编号稳定打破平局

### 第四步：截成 Top-K

代码里的有效保留数是：

```cpp
const int keep_k = std::min(std::max(1, p.CHASE_TOPK_KEEP), NTEST);
```

然后：

- 如果 `good` 候选数大于 `keep_k`
  - 就裁成前 `keep_k` 个
- 如果 `good` 候选数小于等于 `keep_k`
  - 就全部保留

所以最终真正参与后续计算的候选数是：

```text
实际保留数 = min(good候选数, min(max(1, CHASE_TOPK_KEEP), CHASE_NTEST))
```

可以拆开理解成三层限制：

1. `CHASE_TOPK_KEEP` 至少按 1 处理
2. 保留数不会超过 `CHASE_NTEST`
3. 更不会超过当前实际存在的 `good` 候选数

## 4. Top-K 保留之后，后面哪些步骤只看这 K 个

在 `chase_topk_pruned` 里，下面两个步骤都只在保留下来的 `kept_comps` 上做：

- `ML` 选择
- 每个 bit 的 `Cplus / Cminus` 搜索

也就是说，Top-K 裁剪不是只影响最后一小步，而是会同时影响：

- 选哪条码字当 `ML`
- 每个 bit 的外信息 `omega[j]` 怎么算

## 5. 如果 `good` 候选太少，会发生什么

### 情况 A：`good` 候选少于 K

允许。

这时不会强行补满 K，而是：

- 有几个 `good`，就保留几个

### 情况 B：一个 `good` 候选都没有

也允许。

这时：

- `kept_comps` 为空
- `ML` 选不出来

代码会退回到：

- 用信道硬判 `hard_ch`
- 再补一个扩展 parity

并打印一条 warning。

这说明 `chase_topk_pruned` 的 fallback 不是“放弃解码”，而是：

- 当 Top-K good 候选集为空时，退回到 channel hard decision 基线

## 6. 一个具体例子

假设当前参数是：

- `CHASE_NTEST = 64`
- `CHASE_TOPK_KEEP = 8`

并且这次 64 个候选里：

- 只有 13 个是 `good`

那么流程就是：

1. 先生成全部 64 个候选
2. 筛出 13 个 `good`
3. 对这 13 个按 `score` 排序
4. 只保留前 8 个
5. 后续 `ML` 和逐 bit 的 `Cplus/Cminus` 都只在这 8 个里找

再比如：

- `CHASE_NTEST = 64`
- `CHASE_TOPK_KEEP = 64`
- 实际 `good` 只有 11 个

那么最后实际保留的也只是：

- 11 个

不是 64 个。

## 7. 总结

`kChaseTopkKeep` 的选择逻辑可以总结成一句话：

- 在顶层它是 `chase_topk_pruned` 的保留数参数；
- 在 `ofec_sweep` 里候选列表优先于单值；
- 在真正解码时，它不是从“全部候选”里截 Top-K，而是先筛 `good`，再从 `good` 里保留 Top-K；
- 最终实际保留数还会再受 `CHASE_NTEST` 和 `good` 候选总数限制。

所以当前实现下，`kChaseTopkKeep` 更准确地说是：

- `chase_topk_pruned` 中“参与 ML/外信息计算的 good 候选最大保留数”

而不是：

- “从所有候选里无条件取前 K 个”
