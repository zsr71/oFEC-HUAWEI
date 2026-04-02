# MUX 参数之间的关系说明

本文说明下面这几个顶层参数在 `ofec_single` / `ofec_sweep` 中各自控制什么，以及它们之间是如何互相影响的。

```cpp
static const std::vector<int> kSisoActiveList = {64, 64, 64, 64};
static constexpr int  kMuxGroupG          = 1;
static constexpr int  kMuxSchedulingMode  = 0;
static constexpr int  kMuxPriorityRule    = 0;
static constexpr bool kMuxEnableReconfig  = false;
static constexpr int  kMuxBypassScheme    = 1;
```

## 1. 先看整体顺序

在一个软 tile 里，当前调度大致按下面顺序发生：

1. 先做 early-stop 判定
2. 把每个 row 分成：
   - `EarlyStopped`
   - `NeedSiso`
3. 再看当前 tile 的 `SISO` 预算是否充足
4. 如果预算不足，再决定用哪种 MUX 方式去裁 `NeedSiso`
5. 如果预算充足，就基本不需要 MUX 裁剪

所以这几个参数里，真正最先决定“这一轮需不需要调度”的，是：

- `kSisoActiveList`

而不是别的几个。

## 2. `kSisoActiveList` 的作用

`kSisoActiveList[t]` 表示：

- 第 `t` 个 tile 本轮最多允许多少行进入正常 SISO/Chase 解码

它控制的是 **预算大小**。

例如：

```cpp
kSisoActiveList = {64, 64, 64, 64}
```

表示四个 tile 的预算都很大。

如果某个 tile 实际待处理行数不超过这个预算，那么：

- 不需要真正裁剪 `NeedSiso`
- MUX 分组、排序、旁路这些机制即使配置了，也可能几乎看不出效果

所以：

- `kSisoActiveList` 决定“需不需要竞争”
- 后面其他 MUX 参数决定“如果要竞争，该怎么竞争”

## 3. `kMuxGroupG` 的作用

`kMuxGroupG` 控制：

- 预算是按整个 tile 一起分配
- 还是先把 tile 切成若干组后，再按组分配

### 当 `kMuxGroupG = 1`

表示：

- 整个 tile 是一个全局池
- 所有 `NeedSiso` 行放在一起看

这时：

- 如果是旧 MUX，就按原始顺序裁掉后面的
- 如果是新 MUX，就在整个 tile 范围内排序后裁掉后面的

### 当 `kMuxGroupG > 1`

表示：

- 先把当前 tile 的 code 均匀分成 `G` 组
- 再把总 SISO 预算均匀分到每组
- 每组内部各自裁剪

这会直接影响：

- `kMuxSchedulingMode = 0/1` 的排序范围
- `kMuxEnableReconfig = true` 时 staged 调度的分组方式

所以：

- `kMuxGroupG` 决定“在多大范围内竞争”
- 不是决定“用什么排序规则”

## 4. `kMuxSchedulingMode` 的作用

这个参数决定：

- 当需要裁 `NeedSiso` 时，用旧逻辑还是新逻辑

### `kMuxSchedulingMode = 0`

表示使用旧 MUX：

- 不看 early-stop 的细节信息
- 只看当前状态是不是 `NeedSiso`
- 然后按原始顺序保留前面的、裁掉后面的

这是当前的基线调度方式。

### `kMuxSchedulingMode = 1`

表示使用新 MUX：

- 对 `NeedSiso` 行先打分
- 打分只用 early-stop 条件 1 输出的细节信息
- 再按分数从高到低保留

新 MUX 当前只会影响：

- `NeedSiso` 之间谁先留下

不会影响：

- `EarlyStopped` 行

所以：

- `kMuxSchedulingMode` 决定“裁剪时用旧规则还是新规则”
- 但前提是预算真的不够

## 5. `kMuxPriorityRule` 的作用

这个参数只有在：

- `kMuxSchedulingMode = 1`

时才有意义。

如果 `kMuxSchedulingMode = 0`，它基本不生效。

### `kMuxPriorityRule = 0`

表示：

- 更差的 row 优先

当前打分更看重：

- syndrome 里非 0 的个数更多
- overall parity 没过

也就是：

- 越“错得明显”的 row 越优先留下来解

### `kMuxPriorityRule = 1`

表示：

- 更接近通过 early-stop 的 row 优先

当前更偏向：

- BCH 已经通过、只差 overall parity 的 row
- syndrome 非 0 更少的 row

也就是：

- 越“接近过线”的 row 越优先留下来解

所以：

- `kMuxPriorityRule` 不是独立生效参数
- 它是 `kMuxSchedulingMode = 1` 下的子规则

## 6. `kMuxEnableReconfig` 的作用

这个参数决定：

- 是否启用“重配置版 MUX 调度”

### `kMuxEnableReconfig = false`

表示：

- 不走 staged reconfig 调度
- 只做普通预算裁剪

此时真正生效的是：

- `kSisoActiveList`
- `kMuxGroupG`
- `kMuxSchedulingMode`
- `kMuxPriorityRule`（仅当 mode=1）

而：

- `kMuxBypassScheme`

基本不会起实际作用。

### `kMuxEnableReconfig = true`

表示：

- 进入 staged reconfig 调度
- 允许使用额外的 bypass 边来做更复杂的调度

此时：

- `kMuxBypassScheme` 才开始有意义

另外，当前实现里如果：

- `kMuxSchedulingMode = 1`

那么会先按新 MUX 规则对 `NeedSiso` 做一轮预裁剪，
再把留下的结果交给 reconfig 调度继续处理。

所以：

- `kMuxEnableReconfig` 决定“是否启用更复杂的后续调度”
- 并不替代 `kMuxSchedulingMode`

## 7. `kMuxBypassScheme` 的作用

这个参数控制：

- 采用哪一组预定义的 bypass 边集合

例如：

- `1 = scheme1`
- `2 = scheme2`
- 可能还有后续新增的 `scheme3`

但它有一个重要前提：

- 只有 `kMuxEnableReconfig = true` 时，它才真正影响结果

如果：

- `kMuxEnableReconfig = false`

那么虽然这个编号会被传进配置，
但实际调度并不会用到这组 bypass edges。

所以：

- `kMuxBypassScheme` 依赖 `kMuxEnableReconfig`

## 8. 这几个参数的依赖关系

可以把它们的关系简化成下面这样：

### 第一层：预算是否紧张

由：

- `kSisoActiveList`

决定。

如果预算够：

- 后面的 MUX 机制几乎不明显

如果预算不够：

- 才进入真正的调度竞争

### 第二层：竞争范围怎么划

由：

- `kMuxGroupG`

决定。

它控制：

- 全局池竞争
- 还是分组竞争

### 第三层：竞争时按什么规则挑人

由：

- `kMuxSchedulingMode`
- `kMuxPriorityRule`

共同决定。

其中：

- `kMuxPriorityRule` 依赖 `kMuxSchedulingMode = 1`

### 第四层：是否启用重配置和旁路

由：

- `kMuxEnableReconfig`
- `kMuxBypassScheme`

共同决定。

其中：

- `kMuxBypassScheme` 依赖 `kMuxEnableReconfig = true`

## 9. 几组典型配置怎么理解

### 配置 A

```cpp
kSisoActiveList   = {64,64,64,64}
kMuxGroupG        = 1
kMuxSchedulingMode= 0
kMuxPriorityRule  = 0
kMuxEnableReconfig= false
kMuxBypassScheme  = 1
```

含义：

- 预算很宽
- 全局池
- 旧 MUX
- 不启用 reconfig

实际效果通常是：

- 几乎看不到 MUX 的裁剪作用

### 配置 B

```cpp
kSisoActiveList   = {16,16,16,16}
kMuxGroupG        = 1
kMuxSchedulingMode= 1
kMuxPriorityRule  = 0
kMuxEnableReconfig= false
```

含义：

- 预算紧张
- 全局池
- 用新 MUX
- 规则是“更差优先”

实际效果：

- `NeedSiso` 在全 tile 范围内按 syndrome/overall 信息排序

### 配置 C

```cpp
kSisoActiveList   = {16,16,16,16}
kMuxGroupG        = 2
kMuxSchedulingMode= 1
kMuxPriorityRule  = 1
kMuxEnableReconfig= true
kMuxBypassScheme  = 2
```

含义：

- 预算紧张
- 先分两组
- 组内按“更接近通过优先”预裁剪
- 再启用 reconfig + bypass 继续调度

这是最复杂的一种组合。

## 10. 一句话总结

这几个参数之间的关系可以记成：

- `kSisoActiveList` 决定有没有资源竞争
- `kMuxGroupG` 决定竞争范围
- `kMuxSchedulingMode` 决定用旧 MUX 还是新 MUX
- `kMuxPriorityRule` 只在新 MUX 下有效
- `kMuxEnableReconfig` 决定是否启用更复杂的重配置调度
- `kMuxBypassScheme` 只在 reconfig 开启时有效

如果只想先做“旧 MUX vs 新 MUX”的 BER 对比，最简单的做法是先固定：

- `kMuxGroupG = 1`
- `kMuxEnableReconfig = false`

只扫：

- `kMuxSchedulingMode`
- `kMuxPriorityRule`

这样最容易把 MUX 排序策略本身的影响看清楚。
