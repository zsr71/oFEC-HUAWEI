#!/usr/bin/env python3
"""绘制 K>8（GT）调用与未调度 code 的逐调用散点图。"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


DEFAULT_ROUNDS_CSV = Path(
    "data/level56_schedule/ofec_single_level56_schedule_rounds.csv"
)
DEFAULT_CODES_CSV = Path(
    "data/level56_schedule/ofec_single_level56_schedule_codes.csv"
)
DEFAULT_OUTPUT = Path(
    "data/level56_schedule/figures/level56_gt_unscheduled_scatter_3p06dB.png"
)
DEFAULT_FONT = Path(
    "data/level56_schedule/figures/fonts/NotoSansCJKsc-Regular.otf"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="绘制每次五六级共享调用的 GT 状态与未调度 code 数。"
    )
    parser.add_argument("--rounds-csv", type=Path, default=DEFAULT_ROUNDS_CSV)
    parser.add_argument("--codes-csv", type=Path, default=DEFAULT_CODES_CSV)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--font", type=Path, default=DEFAULT_FONT)
    parser.add_argument("--ebn0", type=str, default="3.06 dB")
    return parser.parse_args()


def configure_chinese_font(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(
            f"未找到中文字体：{path}。请通过 --font 指定支持中文的字体。"
        )
    font_manager.fontManager.addfont(str(path))
    family = font_manager.FontProperties(fname=str(path)).get_name()
    plt.rcParams["font.family"] = family
    plt.rcParams["axes.unicode_minus"] = False


def read_invocation_branches(path: Path) -> dict[int, str]:
    branches: dict[int, str] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            invocation = int(row["invocation"])
            branches.setdefault(invocation, row["branch"])
    return branches


def read_unscheduled_counts(path: Path) -> dict[int, int]:
    counts: dict[int, int] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            invocation = int(row["invocation"])
            counts.setdefault(invocation, 0)
            if row["final_action"] == "Unscheduled":
                counts[invocation] += 1
    return counts


def write_invocation_summary(
    path: Path,
    branches: dict[int, str],
    unscheduled: dict[int, int],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["invocation", "branch", "is_gt", "unscheduled_codes"])
        for invocation in sorted(branches):
            branch = branches[invocation]
            writer.writerow(
                [
                    invocation,
                    branch,
                    int(branch == "K_GT_8"),
                    unscheduled.get(invocation, 0),
                ]
            )


def plot_scatter(
    branches: dict[int, str],
    unscheduled: dict[int, int],
    output: Path,
    ebn0: str,
) -> None:
    if not branches:
        raise ValueError("调度 CSV 中没有 invocation 数据")

    gt_x: list[int] = []
    gt_y: list[int] = []
    other_x: list[int] = []
    other_y: list[int] = []
    for invocation in sorted(branches):
        count = unscheduled.get(invocation, 0)
        if branches[invocation] == "K_GT_8":
            gt_x.append(invocation)
            gt_y.append(count)
        else:
            other_x.append(invocation)
            other_y.append(count)

    total_unscheduled = sum(gt_y) + sum(other_y)
    gt_unscheduled = sum(gt_y)
    gt_share = 100.0 * gt_unscheduled / total_unscheduled

    plt.rcParams.update(
        {
            "font.size": 13,
            "axes.titleweight": "bold",
            "axes.labelweight": "bold",
            "axes.edgecolor": "#475569",
            "axes.linewidth": 1.0,
        }
    )

    figure, axis = plt.subplots(figsize=(16, 9), dpi=200)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("#f8fafc")

    axis.scatter(
        other_x,
        other_y,
        s=12,
        color="#64748b",
        alpha=0.35,
        linewidths=0,
        label=f"非 GT 调用（{len(other_x):,} 次）",
        rasterized=True,
        zorder=2,
    )
    axis.scatter(
        gt_x,
        gt_y,
        s=24,
        color="#d33f33",
        alpha=0.78,
        edgecolors="white",
        linewidths=0.25,
        label=f"K>8（GT）调用（{len(gt_x):,} 次）",
        rasterized=True,
        zorder=3,
    )

    axis.set_xlim(-50, max(branches) + 50)
    axis.set_ylim(-0.8, max(max(gt_y), max(other_y)) + 2)
    axis.xaxis.set_major_locator(MultipleLocator(1000))
    axis.xaxis.set_minor_locator(MultipleLocator(500))
    axis.yaxis.set_major_locator(MultipleLocator(5))
    axis.grid(axis="y", color="#cbd5e1", linewidth=0.8, alpha=0.75)
    axis.grid(axis="x", which="major", color="#e2e8f0", linewidth=0.7)
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    axis.set_xlabel("五六级共享调用序号 invocation（时间顺序）", labelpad=12)
    axis.set_ylabel("单次调用的未调度 code 数", labelpad=12)
    axis.set_title(
        f"K>8（GT）调用与未调度 code 的时间分布（Eb/N0 = {ebn0}）",
        fontsize=21,
        pad=20,
    )

    summary = (
        f"GT 调用：{len(gt_x):,} / {len(branches):,} 次（{100 * len(gt_x) / len(branches):.2f}%）\n"
        f"GT 产生未调度：{gt_unscheduled:,} / {total_unscheduled:,} code（{gt_share:.2f}%）"
    )
    axis.text(
        0.012,
        0.972,
        summary,
        transform=axis.transAxes,
        va="top",
        ha="left",
        fontsize=13,
        color="#334155",
        linespacing=1.45,
        bbox={
            "boxstyle": "round,pad=0.45",
            "facecolor": "white",
            "edgecolor": "#cbd5e1",
            "alpha": 0.96,
        },
        zorder=5,
    )
    axis.legend(
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="#cbd5e1",
        framealpha=0.96,
        markerscale=1.4,
    )

    figure.text(
        0.5,
        0.025,
        "每个点代表一次真实调用，未进行分箱、滑动平均或窗口平滑。",
        ha="center",
        color="#64748b",
        fontsize=11,
    )
    figure.tight_layout(rect=(0.035, 0.055, 0.98, 0.97))
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    configure_chinese_font(args.font)
    branches = read_invocation_branches(args.rounds_csv)
    unscheduled = read_unscheduled_counts(args.codes_csv)
    summary_path = args.output.with_name(f"{args.output.stem}_points.csv")
    write_invocation_summary(summary_path, branches, unscheduled)
    plot_scatter(branches, unscheduled, args.output, args.ebn0)
    print(f"已保存散点图：{args.output}")
    print(f"已保存逐调用数据：{summary_path}")


if __name__ == "__main__":
    main()
