# oFEC 早停条件/动作顶层可配改造方案

日期：2026-03-18

## 目标

这份文档说明，如何把当前 oFEC 里的 early-stop 机制改造成“顶层可选择早停判定条件”和“顶层可选择早停后执行动作”的形式，并把相关参数统一暴露到顶层，方便：

1. 在 `ofec_single` 顶层直接切换不同 early-stop 条件和动作；
2. 在 `ofec_sweep` 顶层扫描这些条件、动作及其参数；
3. 直接比较不同组合对 BER 的影响。

下面统一使用 `BER` 这个术语。你前面提到的 `br / b2`，这里都按 `BER` 来理解。

## 当前代码现状

当前代码里，early-stop 实际上有两部分：

1. “判定条件”
2. “判定命中后的动作”

但是这两部分在顶层没有解耦，目前主要是通过一个 `detect_mode` 间接控制。

### 1. 当前顶层可配参数

现在顶层只直接暴露了“是否启用早停”和“检测模式”：

- `include/newcode/params.hpp`
- `include/newcode/ofec_single_runner.hpp`
- `include/newcode/ofec_sweep_runner.hpp`

对应字段包括：

- `ENABLE_EARLY_STOP`
- `EARLY_STOP_DETECT_MODE`
- `EARLY_STOP_V2_LLR_ABS_THRESHOLD`
- `EARLY_STOP_V2_MAX_UNRELIABLE_BITS`

相关接线位置：

- `src/ofec_single/ofec_single_params.cpp`
- `src/ofec_sweep/ofec_sweep_runner.cpp`

### 2. 当前“早停判定条件”

当前 tile 级判定入口在：

- `include/newcode/ofec/earlystop/tile_early_stop_stats.hpp`
- `src/rx/ofec/detail/ofec_tile_impl.ipp`

tile 内部先执行：

```cpp
detect_tile_early_stop(prep.lin_matrix, p, p.EARLY_STOP_DETECT_MODE)
```

目前有两个判据：

#### 条件 1：`detect_tile_early_stop_v1`

位置：

- `src/rx/ofec/earlystop/tile_early_stop_stats.ipp`

规则：

1. 把每一行前 `255` 位做硬判决；
2. 检查 BCH syndrome 是否为 0；
3. 检查 overall parity 是否一致；
4. 每一行都通过时，`all_rows_passed=true`。

#### 条件 2：`detect_tile_early_stop_v2`

位置：

- `src/rx/ofec/earlystop/tile_early_stop_stats2.ipp`

规则：

1. 对每一行统计 `abs(LLR) < threshold` 的 bit 数；
2. 若该行不可靠 bit 数 `<= max_unreliable_bits`，则该行通过；
3. 每一行都通过时，`all_rows_passed=true`。

### 3. 当前“命中早停后的动作”

tile 判定之后，并不会直接“整块不解码”，而是：

1. 先构造 `row_passed_flags`
2. 再转成 `mux_state`
3. 然后在行级决定走哪条动作分支

对应位置：

- `src/rx/ofec/detail/ofec_tile_impl.ipp`
- `src/rx/ofec/mux/mux_state_builder.cpp`
- `src/rx/ofec/ofec_row_decoder_core.cpp`

当前行级分支逻辑是：

- `mux_tag == 1`：走 early-stop 行处理
- `mux_tag == 0`：走正常 Chase
- `mux_tag == 2`：本轮不分配 SISO

当前 early-stop 行处理实际走的是：

- `row_early_stop_process_1(...)`

位置：

- `src/rx/ofec/earlystop/row_early_stop_process_1.ipp`

公式：

```cpp
y2 = (lin - lch) + sign(lin) * beta
```

也就是说，当前代码里的“条件”和“动作”并没有顶层分离：

- 你能选 `detect_mode`
- 但不能单独选“命中后是走 `process_1`、`process_2`、还是别的动作”

## 当前结构存在的问题

### 1. 顶层只能选“条件”，不能单独选“动作”

现在你如果想比较：

- `条件1 + 动作1`
- `条件1 + 动作2`
- `条件2 + 动作1`
- `条件2 + 动作2`

是做不到的。

### 2. 参数语义混在一起

例如当前 `beta` 是全局 Chase / fallback 也在用的参数，但 `row_early_stop_process_1()` 也在用它。

这会带来一个问题：

- 你想扫“早停动作的 beta”
- 但实际改到的是“整个 Chase fallback 的 beta”

这两个变量语义并不相同，不应该强绑在同一个参数上。

这里还需要再明确一层：  
你现在真正想做的是把“early-stop 命中后动作里使用的 beta”和“正常 Chase 解码里使用的 beta”彻底分开。

也就是说，后面应该允许独立配置并独立扫描下面两套参数：

1. `Chase / fallback` 使用的 `beta_list`
2. `early-stop action` 使用的 `beta_list`

这两套 `beta` 不应该继续共用当前同一个 `p.beta / p.beta_list`。

### 3. `v1/v2`、`process_1/process_2` 的命名语义太弱

现在代码里的命名更像“实现编号”，不是“行为语义”：

- `v1 / v2`
- `process_1 / process_2`

这不利于后面做参数扫描，也不利于在 CSV 中清晰记录。

## 总体改造原则

我建议这次改造按下面四个原则做。

### 原则 1：把“判定条件”和“命中动作”彻底解耦

也就是：

- 条件只负责回答：这一行/这个 tile 是否认为“可以早停”
- 动作只负责回答：一旦早停，这一行输出什么

### 原则 2：条件参数和动作参数分开命名

例如：

- `LLR threshold`、`max_unreliable_bits` 属于条件参数
- `beta`、`residual scale`、`hard llr mag` 属于动作参数

不要继续把动作参数挂在全局 Chase 参数上。

### 原则 3：内部可以分组，顶层必须好配、好扫

我建议内部在 `Params` 里做语义分组，但在 `ofec_single::Config` 和 `ofec_sweep::SweepParameterConfig` 层面保持“顶层可直接赋值”的形式，方便你：

1. 手工改单次实验参数；
2. 在 sweep 里直接加候选数组做扫描。

### 原则 4：先保持 MUX 状态机不变

当前 `mux_state` 的三态已经够用：

- `NeedSiso = 0`
- `EarlyStopped = 1`
- `Unscheduled = 2`

第一阶段不建议再引入新的状态码。  
也就是说：

- “条件”负责决定哪些行是 `EarlyStopped`
- “动作”负责决定这些 `EarlyStopped` 的行如何构造输出

这样改动最小，便于先把 BER 扫描跑起来。

## 建议的顶层参数模型

## 1. `Params` 中新增两个核心选择项

建议在 `include/newcode/params.hpp` 中增加：

```cpp
int EARLY_STOP_CONDITION_MODE = 1;
int EARLY_STOP_ACTION_MODE = 1;
```

建议的语义：

- `EARLY_STOP_CONDITION_MODE`
  - `0 = disabled`
  - `1 = BCH+overall parity`
  - `2 = LLR reliability`
  - `3 = reserved for hybrid`

- `EARLY_STOP_ACTION_MODE`
  - `0 = bypass early-stop action`
  - `1 = residual + sign * beta`        // 当前 process_1
  - `2 = residual only / scaled`        // 当前 process_2 语义
  - `3 = hard-decode finish`            // 后续可选

说明：

- `ENABLE_EARLY_STOP` 仍可保留，作为总开关；
- 但 `EARLY_STOP_DETECT_MODE` 应逐步被 `EARLY_STOP_CONDITION_MODE` 替代；
- “条件”和“动作”以后不要再共用一个 mode。

## 2. 条件参数

### 条件 1：BCH + overall parity

这一类条件本身需要的参数不多，但为了统一接口，我建议也显式留出来：

```cpp
bool EARLY_STOP_COND_V1_REQUIRE_BCH = true;
bool EARLY_STOP_COND_V1_REQUIRE_OVERALL = true;
```

这样后面如果你要实验：

- 只看 BCH，不看 overall
- 或者只看硬判合法性，不看扩展 parity

就不需要再改内部代码。

### 条件 2：LLR reliability

建议保留并扩展当前参数：

```cpp
float EARLY_STOP_COND_V2_LLR_ABS_THRESHOLD = 0.5f;
int   EARLY_STOP_COND_V2_MAX_UNRELIABLE_BITS = 8;
bool  EARLY_STOP_COND_V2_INCLUDE_OVERALL = true;
```

其中：

- `LLR_ABS_THRESHOLD`：判定“不可靠 bit”的阈值
- `MAX_UNRELIABLE_BITS`：每行允许的不可靠 bit 个数上限
- `INCLUDE_OVERALL`：统计不可靠 bit 时是否把第 `255` 位 overall parity 算进去

后面如果你要再加：

- 只看前 255 位
- 只看信息位区域
- 只看历史位区域

也可以继续在这个条件组里扩展。

## 3. 动作参数

### 动作 1：`residual + sign * beta`

对应当前 `row_early_stop_process_1()`，但建议以后不要继续复用全局 `p.beta`，而是单独拆成：

```cpp
float EARLY_STOP_ACTION_SIGN_BETA = 0.35f;
std::vector<float> EARLY_STOP_ACTION_SIGN_BETA_LIST;
```

这里我要特别说明：

这个 `beta` 不应该归到“条件 1”或“条件 2”下面，  
它应该归到“动作”下面。

原因很简单：

- 同一个条件 1，可以配动作 1、动作 2、动作 3
- 同一个条件 2，也可以配动作 1、动作 2、动作 3

所以 `beta` 是动作参数，不是条件参数。

但为了满足你现在的实验目标，这里还不够。  
还要继续把这套动作 beta 和正常 Chase 的 beta 分离：

```cpp
float CHASE_BETA = 0.35f;
std::vector<float> CHASE_BETA_LIST;

float EARLY_STOP_ACTION_SIGN_BETA = 0.35f;
std::vector<float> EARLY_STOP_ACTION_SIGN_BETA_LIST;
```

如果不想大改字段名，也至少要做到下面这个语义约束：

1. 旧的 `beta / beta_list` 只服务于正常 Chase / fallback
2. early-stop action 另外引入一套独立的 `beta / beta_list`

这样你后面才能真正比较：

1. 固定 `CHASE_BETA_LIST`，只扫 `EARLY_STOP_ACTION_SIGN_BETA_LIST`
2. 固定 `EARLY_STOP_ACTION_SIGN_BETA_LIST`，只扫 `CHASE_BETA_LIST`
3. 两套 `beta_list` 同时变化，但在 CSV 中分别记录

这也是你当前特别关心的点：  
`v1` 对应的 early-stop 行动作所用的 `beta`，必须和正常 Chase 解码里的 `beta` 分离。

### 动作 2：`residual only / scaled`

对应当前 `row_early_stop_process_2()` 的语义，建议独立出自己的参数：

```cpp
float EARLY_STOP_ACTION_RESIDUAL_DIVISOR = 1.0f;
```

以后可以定义成：

```cpp
y2 = (lin - lch) / EARLY_STOP_ACTION_RESIDUAL_DIVISOR
```

而不是直接复用全局 `ALPHA`。

这样你就可以独立比较：

- 正常 Chase 的 `ALPHA`
- early-stop 动作的 residual scaling

### 动作 3：hard-decode finish

这类动作不是这次第一阶段必须落地，但我建议参数位先预留：

```cpp
float EARLY_STOP_ACTION_HARD_LLR_MAG = 1.0f;
```

如果后面要实现：

- “条件命中后不再跑 Chase，直接走一次硬判收尾”

这个参数就能直接接上。

## 顶层入口怎么改

## 1. `ofec_single::Config`

建议在 [include/newcode/ofec_single_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_single_runner.hpp) 中把当前 early-stop 配置扩成下面这种形式：

```cpp
bool enable_early_stop = true;
int early_stop_condition_mode = 1;
int early_stop_action_mode = 1;

float early_stop_cond_v2_llr_abs_threshold = 0.5f;
int early_stop_cond_v2_max_unreliable_bits = 8;
bool early_stop_cond_v2_include_overall = true;

bool early_stop_cond_v1_require_bch = true;
bool early_stop_cond_v1_require_overall = true;

float chase_beta = 0.35f;
std::vector<float> chase_beta_list;

float early_stop_action_sign_beta = 0.35f;
std::vector<float> early_stop_action_sign_beta_list;
float early_stop_action_residual_divisor = 1.0f;
float early_stop_action_hard_llr_mag = 1.0f;
```

然后在 [src/ofec_single/ofec_single_params.cpp](/home/zsr71/projects/newcode/src/ofec_single/ofec_single_params.cpp) 里统一写回 `Params`。

这样以后在 `apps/ofec_single.cpp` 顶层就可以直接写：

```cpp
constexpr int kEarlyStopConditionMode = 2;
constexpr int kEarlyStopActionMode = 1;
constexpr float kEarlyStopActionSignBeta = 0.25f;
constexpr float kEarlyStopCondV2LlrAbsThreshold = 6.0f;
constexpr int kEarlyStopCondV2MaxUnreliableBits = 4;
```

## 2. `ofec_sweep::SweepParameterConfig`

建议在 [include/newcode/ofec_sweep_runner.hpp](/home/zsr71/projects/newcode/include/newcode/ofec_sweep_runner.hpp) 里新增对应的扫描候选：

```cpp
std::vector<int> early_stop_condition_candidates;
std::vector<int> early_stop_action_candidates;

std::vector<float> chase_beta_candidates;
std::vector<float> early_stop_action_sign_beta_candidates;

std::vector<float> early_stop_cond_v2_llr_abs_threshold_candidates;
std::vector<int> early_stop_cond_v2_max_unreliable_bits_candidates;

std::vector<bool> early_stop_cond_v1_require_bch_candidates;
std::vector<bool> early_stop_cond_v1_require_overall_candidates;
```

其中原则是：

- 不 relevant 的参数，如果候选为空，就用当前标量默认值；
- relevant 的参数，如果候选非空，就进入 scenario 笛卡尔积；
- irrelevant 的参数不参与 scenario 扩展，只写进 CSV 作为记录。

举例：

如果当前扫的是：

- `condition=1`
- `action=1`

那应该参与展开的是：

- `early_stop_condition_candidates`
- `early_stop_action_candidates`
- `chase_beta_candidates`
- `early_stop_action_sign_beta_candidates`
- `early_stop_cond_v1_*`

而 `v2` 的 threshold 和 unreliable bits 应该只记录，不参与展开。

这样能避免 scenario 数量无意义爆炸。

这里我建议明确规定：

1. `chase_beta_candidates` 只影响正常 Chase / fallback
2. `early_stop_action_sign_beta_candidates` 只影响 early-stop action

这样才能支持你真正想做的实验：

- `condition=1, action=1`
- 固定 `CHASE_BETA_LIST`
- 单独扫描 `EARLY_STOP_ACTION_SIGN_BETA_LIST`

## 内部实现怎么改

## 1. 判定层：从“detect_mode”改成“condition mode”

建议把当前：

- `detect_tile_early_stop(...)`

保留接口形式，但内部语义改成：

```cpp
evaluate_early_stop_condition(...)
```

至少逻辑上要清晰成：

- `condition_mode=1` -> BCH + overall parity
- `condition_mode=2` -> unreliable bit count

调用位置仍然在：

- `src/rx/ofec/detail/ofec_tile_impl.ipp`

这里第一阶段不需要动 `TileEarlyStopResult` 的结构，只需要保证：

- `row_passed_flags`
- `rows_passed`
- `all_rows_passed`

依然正常产出。

## 2. 动作层：从 `process_1/process_2` 改成语义分发

建议新增一个统一分发函数，例如：

```cpp
apply_early_stop_action(...)
```

内部再按 `EARLY_STOP_ACTION_MODE` 分发到：

- `action_sign_beta`
- `action_residual_only`
- `action_hard_finish`

这样 [src/rx/ofec/ofec_row_decoder_core.cpp](/home/zsr71/projects/newcode/src/rx/ofec/ofec_row_decoder_core.cpp) 里就不再直接写死：

```cpp
row_early_stop_process_1(...)
```

而是改成：

```cpp
apply_early_stop_action(...)
```

这一步是整个改造的关键。  
只有这一步改掉之后，你才能真正比较“同一个条件配不同动作”的 BER 影响。

## 3. MUX 状态保持不变

当前 [src/rx/ofec/mux/mux_state_builder.cpp](/home/zsr71/projects/newcode/src/rx/ofec/mux/mux_state_builder.cpp) 的三态足够继续使用：

- `NeedSiso`
- `EarlyStopped`
- `Unscheduled`

也就是说：

- 条件层只负责把哪些行置成 `EarlyStopped`
- 动作层只负责这些 `EarlyStopped` 行输出什么

第一阶段不建议引入新的 `StateTag`。

## CSV / 日志怎么改

如果后面要比较 BER，CSV 里必须把“条件”和“动作”以及各自参数都记下来。

我建议在 `ofec_sweep` 结果 CSV 中新增至少这些列：

```text
early_stop_condition_mode
early_stop_action_mode
chase_beta
chase_beta_list
early_stop_action_sign_beta
early_stop_action_sign_beta_list
early_stop_action_residual_divisor
early_stop_cond_v2_llr_abs_threshold
early_stop_cond_v2_max_unreliable_bits
early_stop_cond_v1_require_bch
early_stop_cond_v1_require_overall
tile_early_stop_pct
tile_row_early_stop_pct
```

这样后面你才能直接筛选：

- “同一个 action 下，不同 condition 的 BER”
- “同一个 condition 下，不同 action 的 BER”
- “固定 condition/action，只扫 beta 或 threshold”

## 我建议的落地顺序

### 第 1 步：先把参数模型顶层拆开

先改：

- `Params`
- `ofec_single::Config`
- `ofec_sweep::SweepParameterConfig`
- app 顶层常量

目标是先做到：

- 条件可选
- 动作可选
- 参数可配
- sweep 可记日志

但此时内部逻辑还可以先继续映射到旧实现。

这里第一批就应该先把“两套 beta”拆出来：

1. `beta / beta_list` 继续留给 Chase / fallback
2. 新增 `early_stop_action_sign_beta / early_stop_action_sign_beta_list`

即使第一版内部动作还只是 `process_1`，也要先把参数归属拆干净。

### 第 2 步：把行级动作分发抽象出来

把 `row_early_stop_process_1/2` 改成语义动作分发。

这是最关键的一步，因为它决定后面能不能真正做：

- `condition1 + action2`
- `condition2 + action1`

这种组合实验。

### 第 3 步：扩 sweep scenario builder

在 `ofec_sweep` 里把：

- `condition candidates`
- `action candidates`
- 各自参数 candidates

加入 scenario 构建逻辑。

同时把不相关参数从 scenario 笛卡尔积里剔除，避免扫出来太大。

### 第 4 步：补 CSV 和日志

这一层要保证你最后看结果时，不需要再回头猜：

- 这一行 BER 对应的是哪个条件
- 对应的是哪个动作
- 对应的 beta / threshold / unreliable bits 到底是多少

## 兼容性建议

为了降低一次性改动风险，我建议做一个过渡期：

### 过渡期映射

- `EARLY_STOP_DETECT_MODE` 先保留
- 但内部优先读 `EARLY_STOP_CONDITION_MODE`
- 如果新字段未设置，再回退到旧字段

这样：

1. 旧的 `ofec_single` / `ofec_sweep` 配置还能继续跑；
2. 新的条件/动作框架可以逐步接进来；
3. 不需要一次把所有脚本和常量都改完。

## 验收目标

这次改造完成后，至少应该能直接在顶层比较下面几类组合的 BER：

1. `condition=1, action=1`
2. `condition=1, action=2`
3. `condition=2, action=1`
4. `condition=2, action=2`

并且还能扫描：

1. `action=1` 下的 `beta`
2. `condition=2` 下的 `llr_abs_threshold`
3. `condition=2` 下的 `max_unreliable_bits`

如果 CSV 记录完整，就可以直接拿 sweep 结果做横向比较。

## 一句话结论

我的建议不是“只把 `v1/v2` 再多加几个参数”，而是把 early-stop 彻底拆成两层：

1. `condition`: 什么时候判定可以早停
2. `action`: 判定早停之后怎么生成这一行输出

然后把：

- 条件选择
- 动作选择
- 条件参数
- 动作参数

统一拉到顶层，再接到 `ofec_single` 和 `ofec_sweep`，这样你后面就可以真正系统地比较它们对 BER 的影响，而不是继续把不同变量混在一起扫。
