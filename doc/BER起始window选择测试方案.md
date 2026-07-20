# BER起始Window选择测试方案

## 目标
当前 BER 统计在 [ber.cpp](/home/zsr71/projects/newcode/src/rx/ber/ber.cpp) 中采用固定裁剪：

- 丢掉前 `4` 个 window
- 丢掉后 `1` 个 window

现在需要回答的问题是：

- 前面为什么是 `4`
- 这个值是否足够
- 是否应该改成 `3`、`5`，或者按配置动态确定

测试目标不是先改代码逻辑，而是先通过实验判断：

- 从第几个 window 开始，BER 进入稳定区

## 你提出的方案
你现在的想法是：

1. 保持现有 warmup 不变
2. 不直接只看最终整体 BER
3. 改成按 window 逐段统计 BER
4. 看 BER 在时间方向上的分布
5. 判断从裁掉前多少个 window 开始，后面的 BER 会比较平稳

这个方案是合理的，而且是**最直接、最容易解释**的第一步。

它能直接回答：

- 前几个 window 是否明显更差
- BER 是突然稳定，还是缓慢过渡
- 固定裁掉 `4` 个是否保守、激进，还是刚好

## 推荐先做的第一版方案
### 方案A：逐window BER曲线
这是最推荐先做的。

思路是：

- 保持当前 decoder / warmup / early-stop / MUX 配置不变
- 对解码输出的 info bit，不再只算一个总 BER
- 改成按 window 切块统计

每个 window 统计：

- `window_index`
- `bit_start`
- `bit_end`
- `errs`
- `total`
- `ber`

最后导出 CSV，例如：

```csv
window_index,errs,total,ber
0,120,234432,5.12e-4
1,95,234432,4.05e-4
2,70,234432,2.98e-4
...
```

然后在 Matlab 里画：

- `window_index` 横轴
- `ber` 纵轴

这样你能直接看出：

- BER 是从第几个 window 开始趋稳

### 这个方案的优点
- 最直观
- 改动小
- 结果容易汇报
- 很适合先做探索性分析

### 这个方案的局限
- 如果单个 window 内错误数太少，BER 抖动会比较大
- 高 Eb/N0 下，单 window BER 可能非常 noisy

所以这套方案更适合：

- 中等 BER 区间
- 或者先在相对低一点的 Eb/N0 下看趋势

## 更稳一点的替代方案
### 方案B：滑动窗口平均BER
在方案 A 的基础上，再加一个滑动平均。

例如：

- 每 `1` 个 window 算原始 BER
- 再对连续 `3` 个或 `5` 个 window 做滑动平均

这样会得到两条曲线：

- 原始逐 window BER
- 平滑后的 BER

优点：

- 更容易看出稳定区
- 不容易被单个 window 的随机误差误导

缺点：

- 会弱化边界变化
- 需要明确平均窗口长度

### 方案C：前缀裁剪敏感性测试
这个方案不是看每个 window 的 BER，而是直接测试：

- 如果前面裁掉 `k` 个 window，最终整体 BER 会是多少

例如固定测试：

- `cut_prefix = 0,1,2,3,4,5,6,7`

每个 `cut_prefix` 都重新计算一次：

- `ber_after_cut_prefix_k`

输出表类似：

```csv
cut_prefix_windows,errs,total,ber
0,...
1,...
2,...
3,...
4,...
5,...
```

然后观察：

- 当 `k` 增加到某个值后，BER 变化已经很小

优点：

- 直接回答“前面裁几个更合适”
- 结果和最终 BER 口径最接近

缺点：

- 看不到时域细节
- 只能看到累积后的最终效果

### 方案D：逐window + 前缀裁剪联合测试
这是我更推荐的完整方案，适合作为第二阶段。

先做：

- 方案A，找出 BER 进入稳定的候选区域

再做：

- 方案C，验证裁掉 `k` 个 window 后最终 BER 是否稳定

这样可以同时得到：

- 时域上的直观判断
- 最终整体 BER 口径上的验证

## 我对你当前方案的判断
我认为你的方案是对的，而且应该作为**第一版实验**。

原因是：

- 你现在首先缺的是直观观察
- 先看到 BER 在 window 方向上的分布，才能知道后面该怎么裁
- 不需要一开始就搞复杂统计

但我建议不要只停在“逐 window BER 图”。

更稳的流程应该是：

1. 先画逐 window BER
2. 再做滑动平均辅助观察
3. 最后做前缀裁剪敏感性测试

这样结论会更扎实。

## 推荐实现顺序
### 第一步
新增一个独立 app，例如：

- `apps/ofec_ber_window_probe.cpp`

只做一件事：

- 固定一组解码配置
- 导出每个 window 的 BER

### 第二步
在这个 app 里额外输出滑动平均列，例如：

- `ber_ma3`
- `ber_ma5`

### 第三步
同一个 app 再追加一份前缀裁剪敏感性输出，例如：

- `cut_prefix_windows`
- `ber_after_cut`

## 推荐输出文件
建议至少输出两份 CSV：

### 1. 逐window BER
例如：

- `data/ber_probe/window_ber.csv`

表头建议：

```csv
window_index,errs,total,ber,ber_ma3,ber_ma5
```

### 2. 前缀裁剪敏感性
例如：

- `data/ber_probe/prefix_cut_ber.csv`

表头建议：

```csv
cut_prefix_windows,errs,total,ber
```

## 推荐实验口径
为了避免误判，建议：

- 不只测一个 Eb/N0
- 先选 `2~3` 个代表点

例如：

- 一个中等 BER 点
- 一个较低 BER 点
- 一个接近你目标工作点的点

因为：

- 不同 Eb/N0 下，进入稳定区的速度可能不同

## 判断标准建议
最终可以用下面这种工程化标准来决定前缀裁剪量：

- 当从第 `k` 个 window 开始后，
- 后续若干个 window 的 BER 波动已经落在可接受范围内，
- 且 `cut_prefix = k` 与 `cut_prefix = k+1` 的最终整体 BER 差异已经很小，
- 就认为 `k` 是足够的

更直白一点：

- 先看曲线是否平
- 再看整体 BER 是否已经不敏感

## 一句话结论
你的“逐 window 看 BER 分布”的方案是正确的，应该作为第一步。

更完整的推荐是：

1. 逐 window BER
2. 滑动平均 BER
3. 前缀裁剪敏感性测试

这样才能比较稳地确定：

- 到底应该从第几个 window 开始统计 BER。  
