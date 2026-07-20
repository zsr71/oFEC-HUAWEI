# OFEC 当前 Repo 版与朋友版分类逻辑对比

本文对比当前代码里的两套 hybrid fast classifier：

- `RepoFastClassifier`
- `FriendS1S3Classifier`

对应源码位置：

- [src/rx/ofec/hybrid/hybrid_classifier.ipp](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:133)

本文先分别说明两套逻辑各自怎么判断，再总结它们之间的直接差异。

## 1. RepoFastClassifier 的判断逻辑

代码入口：

- [run_repo_fast_classifier_hard_finish(...)](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:133)

输入是一行 `256-bit` 硬判码字。先计算：

- `S0`：256 位 overall parity syndrome
- `S1`
- `S3`
- `S1^3`

然后按下面顺序判断。

### 1.1 `S0 = 0, S1 = 0, S3 = 0`

- 直接判为 `Clean`
- 说明整体是合法码字

对应代码：

- [hybrid_classifier.ipp:156](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:156)

### 1.2 `S0 = 1, S1 = 0, S3 = 0`

- 判为 `ParityOnly`
- 只翻转 overall parity 位
- 然后做最终合法性校验

对应代码：

- [hybrid_classifier.ipp:159](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:159)

### 1.3 `S1 != 0` 且 `S3 = S1^3`

- 视为主体单错模式
- 通过 `log(S1)` 求出错误位位置
- 翻转该主体位
- 再根据 `S0` 区分：
  - `S0 = 1` -> `OneMain`
  - `S0 = 0` -> `OneMainPlusParity`

对应代码：

- [hybrid_classifier.ipp:163](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:163)

### 1.4 `S0 = 0, S1 != 0, S3 != S1^3`

- 进入两错路径
- 直接调用现有 `BCH t=2` 硬解码器
- 只有在以下条件同时成立时才判为 `TwoMain`
  - decoder 返回成功
  - `corrected_errors == 2`
- 解码成功后回填 255 位主体码字
- 重新计算 overall parity
- 再做最终合法性校验

对应代码：

- [hybrid_classifier.ipp:176](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:176)

### 1.5 其他情况

- 一律 `HardFail`

对应代码：

- [hybrid_classifier.ipp:193](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:193)

## 2. FriendS1S3Classifier 的判断逻辑

代码入口：

- [run_friend_fast_classifier_hard_finish(...)](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:197)

它同样从一行 `256-bit` 硬判码字开始，先计算：

- `S0`
- `S1`
- `S3`

但它的判断主线是先按 `S1/S3` 做分类，再决定后续动作。

### 2.1 `S1 = 0`

如果 `S1 = 0`，继续看 `S3`。

#### 2.1.1 `S1 = 0, S3 != 0`

- 直接判为 `HardFail`
- 朋友版语义里，这表示不可能是 `<= 2` bit error

对应代码：

- [hybrid_classifier.ipp:219](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:219)

#### 2.1.2 `S1 = 0, S3 = 0`

- 先认定 BCH 主体无错
- 再根据 `S0` 区分：
  - `S0 = 1` -> `ParityOnly`
  - `S0 = 0` -> `Clean`

对应代码：

- [hybrid_classifier.ipp:224](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:224)

### 2.2 `S1 != 0` 且 `S3 = S1^3`

- 视为单错模式
- 通过 `log(S1)` 求错误位位置
- 翻转主体错位
- 再根据 `S0` 区分：
  - `S0 = 1` -> `OneMain`
  - `S0 = 0` -> `OneMainPlusParity`

对应代码：

- [hybrid_classifier.ipp:231](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:231)

### 2.3 `S1 != 0` 且 `S3 != S1^3`

这时不直接进 `BCH t=2` 硬解码，而是先算：

- `mu = (S1^3 + S3) / S1^3`
- `tr = Trace(mu)`

对应代码：

- [hybrid_classifier.ipp:247](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:247)

#### 2.3.1 `Trace(mu) != 0`

- 直接判为 `HardFail`
- 当前你定义的语义是：这意味着多错

对应代码：

- [hybrid_classifier.ipp:250](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:250)

#### 2.3.2 `Trace(mu) = 0`

- 当前你定义的语义是：这意味着 2 错
- 于是把该行送入现有 `BCH t=2` 硬解码器
- 只有在以下条件同时成立时才判为 `TwoMain`
  - decoder 返回成功
  - `corrected_errors == 2`
- 解码成功后回填 255 位主体码字
- 重新计算 overall parity
- 再做最终合法性校验

对应代码：

- [hybrid_classifier.ipp:255](/home/zsr71/projects/newcode/src/rx/ofec/hybrid/hybrid_classifier.ipp:255)

## 3. 两套逻辑的直接差异

虽然两边最后都会落到：

- `Clean`
- `ParityOnly`
- `OneMain`
- `OneMainPlusParity`
- `TwoMain`
- `HardFail`

但中间的判定路径不一样。

### 3.1 两错入口条件不同

`RepoFastClassifier` 的两错入口是：

- `S0 = 0`
- `S1 != 0`
- `S3 != S1^3`

也就是说，repo 版先要求 overall parity 必须匹配两错模式，再允许进入 `BCH t=2` 硬解码。

`FriendS1S3Classifier` 的两错入口是：

- `S1 != 0`
- `S3 != S1^3`
- `Trace(mu) = 0`

也就是说，朋友版先用 `Trace(mu)` 把“像 2 错”和“像多错”分开，再决定是否进入 `BCH t=2` 硬解码。

### 3.2 对 `S0` 的使用时机不同

`RepoFastClassifier`：

- `S0` 是顶层硬约束的一部分
- 尤其是在两错路径里，`S0 = 0` 是进入条件

`FriendS1S3Classifier`：

- `S0` 不参与 `mu/trace` 判定
- 主要用于区分：
  - `Clean` vs `ParityOnly`
  - `OneMain` vs `OneMainPlusParity`

换句话说，朋友版的“两错 / 多错”判断核心是 `S1/S3/Trace(mu)`，而不是先看 `S0`。

### 3.3 对 `S1 = 0, S3 != 0` 的解释相同，但出发点不同

两边这类情况最终都失败。

但解释角度不同：

- repo 版：这类输入不落入任何可 hard-finish 的已知模式，所以 `HardFail`
- 朋友版：这是不可能属于 `<=2 bit error` 的模式，所以 `HardFail`

最终动作一致，理论解释不同。

### 3.4 单错路径基本一致

两边在以下判断上是一致的：

- `S1 != 0`
- `S3 = S1^3`

然后：

- 通过 `log(S1)` 找错位
- 翻主体位
- 根据 `S0` 决定是否还要翻 overall parity 位

所以单错路径没有本质差异。

### 3.5 朋友版比 repo 版多了一层 `Trace(mu)` 过滤

这是当前最核心的差异。

repo 版在 `S0 = 0, S1 != 0, S3 != S1^3` 时：

- 直接把该行送进 `BCH t=2` 硬解码器

朋友版在 `S1 != 0, S3 != S1^3` 时：

1. 先算 `Trace(mu)`
2. `Trace(mu) != 0` 直接判多错失败
3. `Trace(mu) = 0` 才送 `BCH t=2` 硬解码器

所以朋友版相当于在两错硬解码前多了一道前置筛选。

### 3.6 两边最终都仍然依赖合法性校验

无论是哪一版，只要走到 hard-finish：

- 都会做 `hard_word_valid_256(...)`
- 要求：
  - 255 位 BCH syndrome 为 0
  - 256 位 overall parity 正确

所以两套逻辑虽然前端分支不同，但最终输出 hard-finish 之前都有同样的出口校验。

## 4. 一句话总结

`RepoFastClassifier` 的核心是：

- 先按 `S0/S1/S3` 直接分 `Clean / ParityOnly / OneMain / TwoMain`
- 两错路径由 `S0 = 0` 约束，再直接送 `BCH t=2` 硬解码

`FriendS1S3Classifier` 的核心是：

- 先按 `S1/S3` 和 `Trace(mu)` 判断是 0 错、1 错、2 错还是多错
- 其中 `Trace(mu) = 0` 被视为 2 错，再送 `BCH t=2` 硬解码
- `Trace(mu) != 0` 直接视为多错失败

所以它们最大的差别不是单错和 0 错，而是：

- **repo 版用 `S0` 决定两错路径入口**
- **朋友版用 `Trace(mu)` 决定两错路径入口**
