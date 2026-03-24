# `new_float_only/apps` 运行说明

本文说明如何用命令行构建并运行 [`new_float_only/apps`](/home/zsr71/projects/newcode/new_float_only/apps) 里的 3 个入口程序：

- `ofec_single_float.cpp`
- `ofec_seed_sweep_float.cpp`
- `ofec_alpha_beta_sweep_float.cpp`

## 1. 先进入工程目录

```bash
cd /home/zsr71/projects/newcode/new_float_only
```

## 2. 配置并编译

首次构建：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

如果之前已经配过 `build/`，后续只需要重新编译：

```bash
cmake --build build -j
```

## 3. 运行单次实验 `ofec_single_float`

对应源码：

- [`ofec_single_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_single_float.cpp)

运行命令：

```bash
./build/ofec_single_float
```

用途：

- 跑一次单独的 float plain 链路实验
- 适合看单组参数下的详细行为

默认特点：

- 会调用 `run_single()`
- 默认开启 `DUMP_WORK_LLR`
- 默认开启 `debug_trace`
- 会输出 Chase trace CSV

常见输出：

- 控制台 BER 信息
- 日志文件：`data/run_<timestamp>_single.log`
- work LLR：`data/llr/work_llr_float.txt`
- Chase CSV：`data/chase_csv/bit1048514/`

## 4. 运行固定参数多 seed 扫描 `ofec_seed_sweep_float`

对应源码：

- [`ofec_seed_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_seed_sweep_float.cpp)

运行命令：

```bash
./build/ofec_seed_sweep_float
```

用途：

- 固定一组 decoder 参数
- 更换多组 seed
- 聚合 pre/post BER

默认特点：

- 走共享的 task scheduler
- 当前示例默认跑 6 个 trial
- 每个 trial 用不同的 `(bitgen_seed, channel_seed)`

常见输出：

- 控制台进度和汇总 BER
- summary CSV：`data/seed_sweep_results_<timestamp>.csv`
- per-trial CSV：`data/seed_sweep_trials_<timestamp>.csv`

## 5. 运行 `alpha/beta` 参数扫描 `ofec_alpha_beta_sweep_float`

对应源码：

- [`ofec_alpha_beta_sweep_float.cpp`](/home/zsr71/projects/newcode/new_float_only/apps/ofec_alpha_beta_sweep_float.cpp)

运行命令：

```bash
./build/ofec_alpha_beta_sweep_float
```

用途：

- 扫描多组 `ALPHA_LIST / beta_list`
- 所有参数模式共享同一批 seed
- 把 `(pattern, seed)` 展开成 task 并行执行
- 最后按 pattern 聚合 BER

默认特点：

- 当前示例会生成一组小网格参数
- 当前示例默认是 `16 patterns × 3 shared seeds = 48 tasks`
- 控制台会按 `post-BER` 从好到差打印 pattern

常见输出：

- 控制台进度
- 控制台 pattern 排名
- summary CSV：`data/ofec_alpha_beta_sweep_float_<timestamp>.csv`

## 6. 只编译某一个目标

如果你只想编译其中一个程序，可以用：

```bash
cmake --build build --target ofec_single_float -j
cmake --build build --target ofec_seed_sweep_float -j
cmake --build build --target ofec_alpha_beta_sweep_float -j
```

## 7. 修改实验参数的方法

这 3 个程序当前都没有命令行参数解析，配置方式是：

1. 直接修改对应 `cpp` 文件顶部的常量
2. 重新编译
3. 再运行可执行文件

例如常改的参数有：

- `Eb/N0`
- `CHASE_L`
- `trial_count`
- `bitgen/channel seed`
- `ALPHA_LIST`
- `beta_list`
- `alpha/beta` 扫描网格

## 8. 一个完整示例

从头到尾跑一次 `alpha/beta` 扫描：

```bash
cd /home/zsr71/projects/newcode/new_float_only
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/ofec_alpha_beta_sweep_float
```

从头到尾跑一次固定参数 seed sweep：

```bash
cd /home/zsr71/projects/newcode/new_float_only
cmake --build build -j
./build/ofec_seed_sweep_float
```

从头到尾跑一次单次实验：

```bash
cd /home/zsr71/projects/newcode/new_float_only
cmake --build build -j
./build/ofec_single_float
```
