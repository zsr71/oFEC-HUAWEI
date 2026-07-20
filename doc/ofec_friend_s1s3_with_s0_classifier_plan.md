# OFEC 第三方案设想：FriendS1S3WithS0Classifier

本文描述一个新的 hybrid fast classifier 设想，先称为：

- `FriendS1S3WithS0Classifier`

它的目标是：

- 保留朋友版“用 `S1/S3/Trace(mu)` 识别主体 2 错”的思路
- 同时把 `S0` 作为扩展 `256-bit` 码字的额外约束加进去
- 避免把 `S0=1` 但 `Trace(mu)=0` 的行也直接送进 `BCH t=2` 硬解码

换句话说，它是一个介于当前两版之间的“第三方案”：

- 比 `RepoFastClassifier` 更强调 `Trace(mu)` 的数学判据
- 比 `FriendS1S3Classifier` 更重视 `S0` 的整体奇偶约束

## 1. 设计目标

当前两套方案在“疑似 2 错”的入口上有明显差异：

- `RepoFastClassifier`
  - 二错入口：`S0=0 && S1!=0 && S3!=S1^3`
- `FriendS1S3Classifier`
  - 二错入口：`S1!=0 && S3!=S1^3 && Trace(mu)=0`

第三方案希望把这两种约束合在一起：

- 主体上必须“像 2 错”
  - 通过 `Trace(mu)=0` 来判断
- 扩展码字整体上也必须“像偶数错”
  - 通过 `S0=0` 来判断

所以第三方案的核心变化只有一条：

- 原朋友版的二错入口：
  - `S1!=0 && S3!=S1^3 && Trace(mu)=0`
- 第三方案的二错入口：
  - `S1!=0 && S3!=S1^3 && Trace(mu)=0 && S0=0`

## 2. 完整判断流程

### 2.1 输入

输入是一行 `256-bit` 硬判码字，先计算：

- `S0`：256 位 overall parity syndrome
- `S1`
- `S3`

若 `S1 != 0`，还需要计算：

- `S1^2`
- `S1^3`

若走到非单错分支，还需要计算：

- `mu = (S1^3 + S3) / S1^3`
- `tr = Trace(mu)`

### 2.2 输出

输出语义仍与当前工程保持一致：

- `Clean`
- `ParityOnly`
- `OneMain`
- `OneMainPlusParity`
- `TwoMain`
- `HardFail`

其中：

- 返回 `true` 表示可以直接 `HardFinish`
- 返回 `false` 表示不能前置硬完成，后续交给 soft path

## 3. 伪代码

```text
function run_friend_s1s3_with_s0_classifier(lin_vec):
    cw = hard_decision_bits_256(lin_vec)

    S0 = overall_parity_syndrome_256(cw)
    [S1, _, S3, _] = bch_255_239_syndromes_1_4_cw_255(cw)

    if S1 == 0:
        if S3 != 0:
            return HardFail, false

        if S0 == 1:
            flip overall parity bit
            if hard_word_valid_256(cw):
                return ParityOnly, true
            else:
                return HardFail, false

        if hard_word_valid_256(cw):
            return Clean, true
        else:
            return HardFail, false

    S1_2 = gf_mul(S1, S1)
    S1_3 = gf_mul(S1_2, S1)

    if S3 == S1_3:
        pos = gf_log(S1)
        if pos invalid:
            return HardFail, false

        flip cw[pos]

        if S0 == 0:
            flip overall parity bit
            if hard_word_valid_256(cw):
                return OneMainPlusParity, true
            else:
                return HardFail, false
        else:
            if hard_word_valid_256(cw):
                return OneMain, true
            else:
                return HardFail, false

    numerator = S1_3 xor S3
    mu = gf_div(numerator, S1_3)
    tr = trace_to_gf2(mu)

    if tr != 0:
        return HardFail, false

    # 关键改动：不仅要求 tr==0，还要求 S0==0
    if S0 != 0:
        return HardFail, false

    ok, decoded, corrected_errors = bch_255_239_decode_hiho_cw_255(cw)
    if not ok:
        return HardFail, false

    if corrected_errors != 2:
        return HardFail, false

    copy decoded[0:255] back to cw
    recompute overall parity bit

    if hard_word_valid_256(cw):
        return TwoMain, true
    else:
        return HardFail, false
```

## 4. 分支解释

### 4.1 `S1 = 0`

这一段完全沿用朋友版当前思路。

- `S1=0, S3!=0`
  - 视为不可能是 `<=2` bit error
  - 直接失败
- `S1=0, S3=0`
  - 说明 BCH 主体无错
  - 再用 `S0` 区分：
    - `S0=0` -> `Clean`
    - `S0=1` -> `ParityOnly`

### 4.2 `S3 = S1^3`

这一段也沿用朋友版当前思路。

- 说明主体符合单错模式
- 通过 `log(S1)` 定位错误位
- 翻主体位
- 再用 `S0` 区分：
  - `S0=1` -> `OneMain`
  - `S0=0` -> `OneMainPlusParity`

### 4.3 `S1!=0 && S3!=S1^3`

这时说明：

- 不是 clean
- 不是 parity-only
- 不是单错

于是进入“可能是 2 错或多错”的分支。

第三方案的判断顺序是：

1. 先算 `Trace(mu)`
2. `tr!=0`
   - 视为多错
   - 失败
3. `tr==0`
   - 先说明“主体上像 2 错”
4. 然后再检查 `S0`
5. `S0!=0`
   - 虽然主体上像 2 错，但整体奇偶不符合偶数错特征
   - 不允许进 `t=2` decoder
   - 失败
6. `S0==0`
   - 这时才真正进入 `BCH t=2` 硬解码

## 5. 与现有两版的关系

### 5.1 相对 RepoFastClassifier

Repo 版的二错入口是：

- `S0=0`
- `S1!=0`
- `S3!=S1^3`

第三方案比它多了一层：

- `Trace(mu)=0`

所以第三方案在二错入口上比 repo 版更严格。

可以理解为：

- repo 版：只要像“非单错且偶数错”，就尝试 `t=2` 解码
- 第三方案：还要再要求 `Trace(mu)=0`，才尝试 `t=2` 解码

### 5.2 相对 FriendS1S3Classifier

朋友版的二错入口是：

- `S1!=0`
- `S3!=S1^3`
- `Trace(mu)=0`

第三方案比它多了一层：

- `S0=0`

所以第三方案在二错入口上比朋友版更保守。

可以理解为：

- 朋友版：只要主体上像 2 错，就尝试 `t=2` 解码
- 第三方案：主体上像 2 错还不够，整体奇偶也必须像偶数错

## 6. 对 compare 结果的预期影响

结合当前 `ofec_fast_classifier_compare.csv` 的结果，第三方案最直接会影响的是这类样本：

- `Trace(mu)=0`
- `corrected_errors=2`
- 但 `S0=1`

例如当前 CSV 里的：

- `lane 0`
- `lane 24`

这些行在当前朋友版下会被当成 `TwoCandidate`，并进一步送入 `BCH t=2` 解码；
在第三方案下，这些行会在进入 `t=2` 解码之前先被 `S0!=0` 拦下来。

也就是说，第三方案的直观效果是：

- 保留朋友版对大部分 `S0=0, tr=0` 的二错识别能力
- 去掉 `S0=1, tr=0` 这批边界样本的前置硬分流

## 7. 这个方案的优缺点

### 优点

- 保留了朋友版基于 `Trace(mu)` 的二错数学判据
- 用上了扩展码字的 `S0` 信息
- 对“疑似 2 错”的放行更稳健
- 更容易解释为“主体判据 + 整体奇偶约束”双重门控

### 风险

- 可能会错过一部分 `Trace(mu)=0` 且 `BCH t=2` 实际可解、但 `S0!=0` 的样本
- 也就是说，它比朋友版更保守，硬分流收益可能会下降

## 8. 一句话总结

第三方案的核心就是：

- 零错和单错分支保持朋友版不变
- 只在二错入口上新增 `S0=0` 约束

即：

- **朋友版：`Trace(mu)=0` 就进两错硬解码**
- **第三方案：`Trace(mu)=0` 且 `S0=0` 才进两错硬解码**
