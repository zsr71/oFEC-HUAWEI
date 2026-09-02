# TwoMain 实验框架

本目录保存 TwoMain 候选方案的受版本控制实验定义、基线签名和小型结论。大型日志、逐码字调度表和完整错误位置保存在仓库根目录下的 `runs/twomain/`，不进入普通 Git 历史。

## 一、当前阶段的边界

阶段二建立实验框架和冻结基线。阶段三已增加方案六统一参数化输出能力并完成固定长帧单随机种子筛选，但默认仍为Legacy；分类、调度、BCH纠正、写回坐标和history语义没有改变。当前可运行的三组基线是：

| 编号 | 含义 | HISO 类别掩码 | 五六级资源 |
|---|---|---:|---:|
| `TM-A0` | 原始方案，所有现有类别均允许使用 HISO | `15` | `8 HISO + 8 SISO` |
| `TM-A1` | 只禁止 TwoMain 使用 HISO | `7` | `8 HISO + 8 SISO` |
| `TM-A2` | 纯 SISO 参考，所有类别均禁止使用 HISO | `0` | `0 HISO + 8 SISO` |

三组正式配置固定使用同一长帧、同一信道和同一 FIFO 条件，以保证结果可直接比较。历史签名分别为 `159`、`9`、`11` 个 Post-FEC 错误，比较区域均为 `58,608,000` bit，且 `forced_evicted=0`。

## 二、目录结构

```text
experiments/twomain/
├── README.md
├── registry.md
├── run_manifest_fields.md
├── specs/
│   ├── baseline_legacy.toml
│   ├── baseline_no_twomain_hiso.toml
│   ├── baseline_pure_siso.toml
│   ├── regression_default_quick.toml
│   ├── scheme6_uniform_32_smoke.toml
│   └── scheme6_differential_32_smoke.toml
└── summaries/
    └── README.md
```

运行后生成：

```text
runs/twomain/<experiment_id>/<run_id>/
├── effective_config.toml
├── manifest.json
├── stdout.log
├── result.json
├── post_fec_error_positions.txt
└── data/
    ├── level56_schedule/
    ├── early_stop_hist/
    └── early_stop_debug/
```

## 三、运行方法

先构建单点程序：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target ofec_single -j
```

查看配置是否能被正确解析，不实际运行：

```bash
python3 scripts/run_twomain_experiment.py \
  --spec experiments/twomain/specs/baseline_legacy.toml \
  --dry-run
```

运行正式 Legacy 基线并核对冻结签名：

```bash
python3 scripts/run_twomain_experiment.py \
  --spec experiments/twomain/specs/baseline_legacy.toml \
  --verify-expected
```

运行另外两组正式基线：

```bash
python3 scripts/run_twomain_experiment.py \
  --spec experiments/twomain/specs/baseline_no_twomain_hiso.toml \
  --verify-expected
```

```bash
python3 scripts/run_twomain_experiment.py \
  --spec experiments/twomain/specs/baseline_pure_siso.toml \
  --verify-expected
```

多随机种子实验可复用同一算法配置，并成对覆盖信息比特和信道随机种子。两个参数必须同时给出，替换后的完整配置会记录到运行清单和`effective_config.toml`：

```bash
python3 scripts/run_twomain_experiment.py \
  --spec experiments/twomain/specs/scheme6a_m2_48_long.toml \
  --bitgen-seed 731941 \
  --channel-seed 810271 \
  --run-id stage3_multiseed_seed02_20260901 \
  --allow-dirty
```

快速执行默认行为等价回归：

```bash
python3 scripts/check_twomain_default_equivalence.py
```

正式结果默认要求工作区干净。如果只是开发实验框架或定位问题，可以显式增加 `--allow-dirty`；这种运行会在清单中记录为脏工作区，不能直接用作最终性能结论。

## 四、配置与程序的关系

配置文件记录完整实验条件。其中能由当前 `ofec_single` 在运行时覆盖的字段，由运行器转换为环境变量；其余固定字段用于核对当前代码默认值和保存完整实验语义。

方案六参数已经进入公共解码参数和统一运行入口。配置中的 `output_policy=legacy_fixed` 显式选择原输出路径；`output_policy=scheme6_parameterized` 时，运行器进一步传递 `m2`、`rho_corr` 和 `rho_keep`。未显式启用新策略时，快速回归签名保持不变。

两份方案六短帧配置只用于机制冒烟，不代表正式入选参数。阶段三前四步的定点证据和冒烟结果见：

```text
experiments/twomain/summaries/stage3_first_four_steps_20260901.md
```

方案六-A/B/C的固定长帧单随机种子筛选结果见：

```text
experiments/twomain/summaries/stage3_single_seed_screening_20260901.md
```

六-A五组配对随机种子验证选择`M2=48`为主候选、`M2=64`为相邻幅度对照。六-B方向一、方向二以及旧方向下的六-C均已完成第一轮筛选，均不推进。

五组配对随机种子验证现已完成，结果见：

```text
experiments/twomain/summaries/stage3_multiseed_validation_20260902.md
```

六-A `M2=48`在五个seed中均改善Legacy，累计错误从1061降到688，但仍高于禁止TwoMain HISO的105；因此保留为参考候选，但不能视为已经解决TwoMain问题。`M2=64`不再推进。

六-B现在区分为：

```text
方向一：rho_corr=1，降低rho_keep；已测试，暂停
方向二：rho_keep=1，降低rho_corr；已测试，未入选
```

方向二固定`M2=48`，在seed01上比较`rho_corr=0.75/0.5/0.25`，分别得到85/74/75错；三组均无forced eviction，但均未优于六-A的64错，因此不进入五组随机种子。详细记录见：

```text
experiments/twomain/summaries/stage3_scheme6b_direction2_seed01_20260902.md
```

## 五、结果使用原则

正式结论只引用满足以下条件的运行：

1. 工作区干净；
2. Git 提交明确；
3. 使用受版本控制的配置；
4. 程序正常完成；
5. 配置中的冻结签名核对通过；
6. 没有发生非预期的强制顶出。

快速回归只负责发现默认行为是否被无意改变，不替代三组完整长帧基线，也不用于比较最终 BER。
