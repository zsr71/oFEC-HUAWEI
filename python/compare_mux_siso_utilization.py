#!/usr/bin/env python3
"""Compare SISO utilization between grouped-budget and reconfig schemes.

Default comparison:
- scheme_static:  G=2, reconfig=false
- scheme_reconfig: G=4, reconfig=true

The simulation treats each code/row as an active decode request with some
probability. Early-stopped rows are simply the inactive rows.
"""

from __future__ import annotations

import argparse
import csv
import random
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

from group_scheduler import schedule_scheme_b_reconfig, schedule_scheme_c_staged


@dataclass
class SchemeStats:
    avg_active: float = 0.0
    avg_scheduled: float = 0.0
    avg_waiting: float = 0.0
    avg_utilization: float = 0.0


def build_groups_even(total: int, group_g: int) -> List[Tuple[int, int]]:
    if group_g < 1:
        raise ValueError("group_g must be >= 1")
    if total > 0 and group_g > total:
        raise ValueError("group_g must be <= total when total > 0")

    base = total // group_g
    rem = total % group_g
    groups: List[Tuple[int, int]] = []
    start = 0
    for i in range(group_g):
        span = base + (1 if i < rem else 0)
        groups.append((start, start + span))
        start += span
    return groups


def split_budget_even(total_budget: int, group_g: int) -> List[int]:
    if total_budget < 0:
        raise ValueError("total_budget must be >= 0")
    if group_g < 1:
        raise ValueError("group_g must be >= 1")

    base = total_budget // group_g
    rem = total_budget % group_g
    budgets = [base] * group_g
    for i in range(rem):
        budgets[i] += 1
    return budgets


def simulate_grouped_budget(
    active_codes: Sequence[int],
    n_code: int,
    n_siso: int,
    group_g: int,
) -> Tuple[int, int]:
    """Mirror C++ apply_siso_budget_grouped semantics.

    Returns:
    - scheduled_count
    - waiting_count
    """
    if group_g <= 1:
        scheduled = min(len(active_codes), n_siso)
        return scheduled, len(active_codes) - scheduled

    groups = build_groups_even(n_code, group_g)
    budgets = split_budget_even(n_siso, group_g)
    active_set = set(active_codes)

    scheduled = 0
    for group_idx, (begin, end) in enumerate(groups):
        local_active = [idx for idx in range(begin, end) if idx in active_set]
        scheduled += min(len(local_active), budgets[group_idx])
    waiting = len(active_codes) - scheduled
    return scheduled, waiting


def simulate_reconfig(
    active_codes: Sequence[int],
    n_code: int,
    n_siso: int,
    group_g: int,
) -> Tuple[int, int, int]:
    """Mirror C++ schedule_scheme_b_reconfig_cpp semantics.

    Returns:
    - stage1_scheduled
    - final_scheduled
    - waiting_count
    """
    free_siso = list(range(n_siso))
    stage1_match, final_match, waiting_codes = schedule_scheme_b_reconfig(
        N_code=n_code,
        N_siso=n_siso,
        G=group_g,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    return len(stage1_match), len(final_match), len(waiting_codes)


def simulate_scheme_c(
    active_codes: Sequence[int],
    n_code: int,
    n_siso: int,
    group_g: int,
) -> Tuple[int, int]:
    """Simulate scheme C staged scheduling."""
    free_siso = list(range(n_siso))
    final_match, waiting_codes = schedule_scheme_c_staged(
        N_code=n_code,
        N_siso=n_siso,
        G=group_g,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    return len(final_match), len(waiting_codes)


def random_active_codes(n_code: int, active_prob: float, rng: random.Random) -> List[int]:
    return [idx for idx in range(n_code) if rng.random() < active_prob]


def run_trials(
    n_code: int,
    n_siso: int,
    active_prob: float,
    trials: int,
    seed: int,
    static_g_ref: int,
    static_g_same_as_reconfig: int,
    reconfig_g: int,
) -> Tuple[SchemeStats, SchemeStats, SchemeStats, SchemeStats, float]:
    rng = random.Random(seed)
    static_ref_stats = SchemeStats()
    static_sameg_stats = SchemeStats()
    reconfig_stats = SchemeStats()
    scheme_c_stats = SchemeStats()
    avg_stage1_util = 0.0

    for _ in range(trials):
        active_codes = random_active_codes(n_code, active_prob, rng)
        active_count = len(active_codes)

        scheduled_static_ref, waiting_static_ref = simulate_grouped_budget(
            active_codes, n_code, n_siso, static_g_ref
        )
        scheduled_static_sameg, waiting_static_sameg = simulate_grouped_budget(
            active_codes, n_code, n_siso, static_g_same_as_reconfig
        )
        stage1_reconfig, final_reconfig, waiting_reconfig = simulate_reconfig(
            active_codes, n_code, n_siso, reconfig_g
        )
        scheduled_scheme_c, waiting_scheme_c = simulate_scheme_c(
            active_codes, n_code, n_siso, reconfig_g
        )

        static_ref_stats.avg_active += active_count
        static_ref_stats.avg_scheduled += scheduled_static_ref
        static_ref_stats.avg_waiting += waiting_static_ref
        static_ref_stats.avg_utilization += (
            scheduled_static_ref / n_siso if n_siso > 0 else 0.0
        )

        static_sameg_stats.avg_active += active_count
        static_sameg_stats.avg_scheduled += scheduled_static_sameg
        static_sameg_stats.avg_waiting += waiting_static_sameg
        static_sameg_stats.avg_utilization += (
            scheduled_static_sameg / n_siso if n_siso > 0 else 0.0
        )

        reconfig_stats.avg_active += active_count
        reconfig_stats.avg_scheduled += final_reconfig
        reconfig_stats.avg_waiting += waiting_reconfig
        reconfig_stats.avg_utilization += final_reconfig / n_siso if n_siso > 0 else 0.0

        scheme_c_stats.avg_active += active_count
        scheme_c_stats.avg_scheduled += scheduled_scheme_c
        scheme_c_stats.avg_waiting += waiting_scheme_c
        scheme_c_stats.avg_utilization += scheduled_scheme_c / n_siso if n_siso > 0 else 0.0

        avg_stage1_util += stage1_reconfig / n_siso if n_siso > 0 else 0.0

    if trials > 0:
        for stats in (static_ref_stats, static_sameg_stats, reconfig_stats, scheme_c_stats):
            stats.avg_active /= trials
            stats.avg_scheduled /= trials
            stats.avg_waiting /= trials
            stats.avg_utilization /= trials
        avg_stage1_util /= trials

    return static_ref_stats, static_sameg_stats, reconfig_stats, scheme_c_stats, avg_stage1_util


def parse_probs(raw: str) -> List[float]:
    values = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if not values:
        raise ValueError("at least one active probability is required")
    for value in values:
        if value < 0.0 or value > 1.0:
            raise ValueError("active probabilities must be in [0, 1]")
    return values


def print_table(
    probs: Iterable[float],
    rows: Sequence[Tuple[SchemeStats, SchemeStats, SchemeStats, SchemeStats, float]],
    static_g_ref: int,
    static_g_same_as_reconfig: int,
    reconfig_g: int,
) -> None:
    header = (
        f"激活概率 | 平均活跃code数 | 静态方案利用率(G={static_g_ref},关) "
        f"| 静态方案利用率(G={static_g_same_as_reconfig},关) "
        f"| 重排方案阶段1利用率(G={reconfig_g},开) "
        f"| 重排方案最终利用率(G={reconfig_g},开) "
        f"| 方案C最终利用率(G={reconfig_g},开) "
        f"| 静态方案平均已调度(G={static_g_ref}) "
        f"| 静态方案平均已调度(G={static_g_same_as_reconfig}) "
        f"| 重排方案平均已调度(G={reconfig_g}) "
        f"| 方案C平均已调度(G={reconfig_g}) "
        f"| 静态方案平均未调度(G={static_g_ref}) "
        f"| 静态方案平均未调度(G={static_g_same_as_reconfig}) "
        f"| 重排方案平均未调度(G={reconfig_g}) "
        f"| 方案C平均未调度(G={reconfig_g})"
    )
    print(header)
    print("-" * len(header))
    for active_prob, (static_ref_stats, static_sameg_stats, reconfig_stats, scheme_c_stats, stage1_util) in zip(probs, rows):
        print(
            f"{active_prob:8.3f} | "
            f"{static_ref_stats.avg_active:10.3f} | "
            f"{static_ref_stats.avg_utilization:21.4f} | "
            f"{static_sameg_stats.avg_utilization:21.4f} | "
            f"{stage1_util:24.4f} | "
            f"{reconfig_stats.avg_utilization:24.4f} | "
            f"{scheme_c_stats.avg_utilization:19.4f} | "
            f"{static_ref_stats.avg_scheduled:23.3f} | "
            f"{static_sameg_stats.avg_scheduled:23.3f} | "
            f"{reconfig_stats.avg_scheduled:20.3f} | "
            f"{scheme_c_stats.avg_scheduled:17.3f} | "
            f"{static_ref_stats.avg_waiting:23.3f} | "
            f"{static_sameg_stats.avg_waiting:23.3f} | "
            f"{reconfig_stats.avg_waiting:20.3f} | "
            f"{scheme_c_stats.avg_waiting:17.3f}"
        )


def write_csv(
    output_path: Path,
    probs: Iterable[float],
    rows: Sequence[Tuple[SchemeStats, SchemeStats, SchemeStats, SchemeStats, float]],
    static_g_ref: int,
    static_g_same_as_reconfig: int,
    reconfig_g: int,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "激活概率",
                "平均活跃code数",
                f"静态方案利用率(G={static_g_ref},关)",
                f"静态方案利用率(G={static_g_same_as_reconfig},关)",
                f"重排方案阶段1利用率(G={reconfig_g},开)",
                f"重排方案最终利用率(G={reconfig_g},开)",
                f"方案C最终利用率(G={reconfig_g},开)",
                f"静态方案平均已调度(G={static_g_ref})",
                f"静态方案平均已调度(G={static_g_same_as_reconfig})",
                f"重排方案平均已调度(G={reconfig_g})",
                f"方案C平均已调度(G={reconfig_g})",
                f"静态方案平均未调度(G={static_g_ref})",
                f"静态方案平均未调度(G={static_g_same_as_reconfig})",
                f"重排方案平均未调度(G={reconfig_g})",
                f"方案C平均未调度(G={reconfig_g})",
            ]
        )
        for active_prob, (
            static_ref_stats,
            static_sameg_stats,
            reconfig_stats,
            scheme_c_stats,
            stage1_util,
        ) in zip(probs, rows):
            writer.writerow(
                [
                    f"{active_prob:.6f}",
                    f"{static_ref_stats.avg_active:.6f}",
                    f"{static_ref_stats.avg_utilization:.6f}",
                    f"{static_sameg_stats.avg_utilization:.6f}",
                    f"{stage1_util:.6f}",
                    f"{reconfig_stats.avg_utilization:.6f}",
                    f"{scheme_c_stats.avg_utilization:.6f}",
                    f"{static_ref_stats.avg_scheduled:.6f}",
                    f"{static_sameg_stats.avg_scheduled:.6f}",
                    f"{reconfig_stats.avg_scheduled:.6f}",
                    f"{scheme_c_stats.avg_scheduled:.6f}",
                    f"{static_ref_stats.avg_waiting:.6f}",
                    f"{static_sameg_stats.avg_waiting:.6f}",
                    f"{reconfig_stats.avg_waiting:.6f}",
                    f"{scheme_c_stats.avg_waiting:.6f}",
                ]
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-code", type=int, default=32)
    parser.add_argument("--n-siso", type=int, default=16)
    parser.add_argument("--static-g", type=int, default=2,
                        help="baseline G for grouped budget (reconfig=false)")
    parser.add_argument("--reconfig-g", type=int, default=4,
                        help="G for scheme B reconfig")
    parser.add_argument("--trials", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20260311)
    parser.add_argument(
        "--active-probs",
        type=parse_probs,
        default=parse_probs("0.20,0.35,0.50,0.65,0.80,0.95"),
        help="comma-separated probabilities that a code needs SISO",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("python/compare_mux_siso_utilization_results.csv"),
        help="path to output CSV table",
    )
    args = parser.parse_args()

    if args.n_code <= 0 or args.n_siso < 0:
        raise ValueError("n_code must be > 0 and n_siso must be >= 0")
    if args.trials <= 0:
        raise ValueError("trials must be > 0")

    print(
        f"对比方案: 静态分组(G={args.static_g}, reconfig=false) vs "
        f"两阶段重排(G={args.reconfig_g}, reconfig=true)"
    )
    print(
        f"附加对照: 静态分组(G={args.reconfig_g}, reconfig=false)"
    )
    print(
        f"新增对照: 方案C顺序式两阶段调度(G={args.reconfig_g}, reconfig=true)"
    )
    print(
        f"参数: n_code={args.n_code}, n_siso={args.n_siso}, "
        f"trials={args.trials}, seed={args.seed}"
    )
    print()

    rows = [
        run_trials(
            n_code=args.n_code,
            n_siso=args.n_siso,
            active_prob=active_prob,
            trials=args.trials,
            seed=args.seed + idx,
            static_g_ref=args.static_g,
            static_g_same_as_reconfig=args.reconfig_g,
            reconfig_g=args.reconfig_g,
        )
        for idx, active_prob in enumerate(args.active_probs)
    ]
    write_csv(args.output, args.active_probs, rows, args.static_g, args.reconfig_g, args.reconfig_g)
    print_table(args.active_probs, rows, args.static_g, args.reconfig_g, args.reconfig_g)
    print()
    print(f"结果表已写入: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
