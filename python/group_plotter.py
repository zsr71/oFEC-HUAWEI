#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_plotter.py
4 组拓扑下方案 A / B 的可视化模块。
这里只负责画图，不负责调度计算。
"""

from typing import Dict, List, Sequence, Tuple

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


def configure_matplotlib_fonts() -> None:
    """优先选择常见中文字体，避免图中中文注释显示为方框。"""
    candidate_fonts = [
        "Microsoft YaHei",
        "SimHei",
        "Noto Sans CJK SC",
        "Source Han Sans SC",
        "WenQuanYi Zen Hei",
        "Arial Unicode MS",
    ]
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}

    for font_name in candidate_fonts:
        if font_name in available_fonts:
            plt.rcParams["font.sans-serif"] = [font_name] + list(plt.rcParams.get("font.sans-serif", []))
            break

    plt.rcParams["axes.unicode_minus"] = False


def y_pos(i: int, total: int) -> float:
    """把离散节点编号映射到 [0, 1] 纵轴坐标，编号越小越靠上。"""
    return 1.0 - (i + 0.5) / total


def group_rect_params(start: int, size: int, total: int) -> Tuple[float, float]:
    """返回覆盖一整组节点的矩形参数 (y0, h)。"""
    end = start + size
    y_top = 1.0 - start / total
    y_bottom = 1.0 - end / total
    return y_bottom, (y_top - y_bottom)


def chunk_labels(prefix: str, indices: Sequence[int], chunk_size: int = 8) -> List[str]:
    """把较长的编号列表拆成多行文本，便于放进图中说明区域。"""
    labels = [f"{prefix}{idx + 1}" for idx in indices]
    if not labels:
        return ["None"]
    return [", ".join(labels[i:i + chunk_size]) for i in range(0, len(labels), chunk_size)]


def format_matching_lines(code_to_siso: Dict[int, int], chunk_size: int = 4) -> List[str]:
    """把匹配结果格式化成多行，便于图中展示。"""
    pairs = [f"C{code_idx + 1}-S{siso_idx + 1}" for code_idx, siso_idx in sorted(code_to_siso.items())]
    if not pairs:
        return ["None"]
    return [", ".join(pairs[i:i + chunk_size]) for i in range(0, len(pairs), chunk_size)]


def draw_scheme_a_schedule(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
    edges: Sequence[Tuple[int, int]],
    match: Dict[int, int],
    waiting_codes: Sequence[int],
    save: str | None = None,
    show: bool = True,
) -> None:
    """绘制方案 A：全局最大匹配。"""
    configure_matplotlib_fonts()

    active_codes = sorted(set(active_codes))
    free_siso = sorted(set(free_siso))
    active_code_set = set(active_codes)
    matched_codes = set(match.keys())
    matched_siso = set(match.values())

    base_edge_color = "#b8b8b8"
    matched_edge_color = "#4f4f4f"

    code_inactive_face = "white"
    code_active_wait_face = "#f2c879"
    code_active_match_face = "#74c69d"
    siso_idle_face = "white"
    siso_busy_face = "#6ea8d9"

    x_code = 0.10
    x_siso = 0.90
    rect_w = 0.14
    rect_pad = 0.02

    fig, ax = plt.subplots(figsize=(12, 9))
    fig.subplots_adjust(left=0.06, right=0.94, top=0.88, bottom=0.23)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    code_xy = [(x_code, y_pos(i, N_code)) for i in range(N_code)]
    siso_xy = [(x_siso, y_pos(j, N_siso)) for j in range(N_siso)]

    code_per_g = N_code // G
    siso_per_g = N_siso // G

    for g in range(G):
        c_start = g * code_per_g
        y0, h = group_rect_params(c_start, code_per_g, N_code)
        ax.add_patch(Rectangle((x_code - rect_w - rect_pad, y0), rect_w, h, fill=False, linewidth=2))

        s_start = g * siso_per_g
        y0s, hs = group_rect_params(s_start, siso_per_g, N_siso)
        ax.add_patch(Rectangle((x_siso + rect_pad, y0s), rect_w, hs, fill=False, linewidth=2))

        ax.text(x_code - rect_w - rect_pad - 0.03, y0 + h / 2, f"G{g}", va="center", ha="right", fontsize=11)
        ax.text(x_siso + rect_pad + rect_w + 0.03, y0s + hs / 2, f"G{g}", va="center", ha="left", fontsize=11)

    for code_idx, siso_idx in edges:
        x1, y1 = code_xy[code_idx]
        x2, y2 = siso_xy[siso_idx]
        ax.plot([x1, x2], [y1, y2], linewidth=0.8, alpha=0.28, color=base_edge_color, zorder=1)

    for code_idx, siso_idx in sorted(match.items()):
        x1, y1 = code_xy[code_idx]
        x2, y2 = siso_xy[siso_idx]
        ax.plot([x1, x2], [y1, y2], linewidth=2.6, alpha=0.95, color=matched_edge_color, zorder=4)

    for code_idx, (x, y) in enumerate(code_xy):
        if code_idx in matched_codes:
            face_color = code_active_match_face
            edge_width = 1.8
        elif code_idx in active_code_set:
            face_color = code_active_wait_face
            edge_width = 1.6
        else:
            face_color = code_inactive_face
            edge_width = 1.2

        ax.scatter(x, y, s=60, facecolors=face_color, edgecolors="black", linewidths=edge_width, zorder=5)
        font_weight = "bold" if code_idx in active_code_set else "normal"
        ax.text(x - 0.02, y, f"C{code_idx + 1}", va="center", ha="right", fontsize=9, fontweight=font_weight)

    free_siso_set = set(free_siso)
    for siso_idx, (x, y) in enumerate(siso_xy):
        if siso_idx in matched_siso:
            face_color = siso_busy_face
            edge_width = 1.8
        elif siso_idx in free_siso_set:
            face_color = siso_idle_face
            edge_width = 1.4
        else:
            face_color = "#dddddd"
            edge_width = 1.2

        ax.scatter(x, y, s=70, facecolors=face_color, edgecolors="black", linewidths=edge_width, zorder=5)
        ax.text(x + 0.02, y, f"S{siso_idx + 1}", va="center", ha="left", fontsize=9)

    max_possible = min(len(active_codes), len(free_siso))
    ax.text(
        0.5,
        1.03,
        (
            "Scheme A: Global Maximum Matching"
            f" | Active Codes={len(active_codes)}"
            f" | Free SISOs={len(free_siso)}"
            f" | Scheduled={len(match)}/{max_possible}"
        ),
        ha="center",
        va="bottom",
        transform=ax.transAxes,
        fontsize=13,
    )

    legend_items = [
        Line2D([0], [0], color=base_edge_color, lw=1.2, label="Allowed edges (topology background)"),
        Line2D([0], [0], color=matched_edge_color, lw=2.8, label="Matched edges in this round"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=code_active_match_face, markersize=8, lw=0, label="Active and matched code"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=code_active_wait_face, markersize=8, lw=0, label="Active but unmatched code"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=siso_busy_face, markersize=8, lw=0, label="Occupied SISO"),
    ]
    ax.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False, fontsize=10)

    active_lines = chunk_labels("C", active_codes, chunk_size=8)
    free_siso_lines = chunk_labels("S", free_siso, chunk_size=8)
    waiting_lines = chunk_labels("C", waiting_codes, chunk_size=8)
    match_lines = format_matching_lines(match, chunk_size=4)

    info_lines = [
        "Scheduling input:",
        f"  Active codes: {active_lines[0]}",
    ]
    info_lines.extend([f"             {line}" for line in active_lines[1:]])
    info_lines.append(f"  Free SISOs: {free_siso_lines[0]}")
    info_lines.extend([f"             {line}" for line in free_siso_lines[1:]])
    info_lines.append(f"  Match count: {len(match)}")
    info_lines.append(f"  Match result: {match_lines[0]}")
    info_lines.extend([f"           {line}" for line in match_lines[1:]])
    info_lines.append(f"  Waiting queue: {waiting_lines[0]}")
    info_lines.extend([f"           {line}" for line in waiting_lines[1:]])

    fig.text(0.07, 0.04, "\n".join(info_lines), ha="left", va="bottom", fontsize=10)

    if save:
        plt.savefig(save, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def draw_scheme_b_schedule(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
    edges: Sequence[Tuple[int, int]],
    stage1_match: Dict[int, int],
    final_match: Dict[int, int],
    waiting_codes: Sequence[int],
    save: str | None = None,
    show: bool = True,
) -> None:
    """绘制方案 B：两阶段 + 可重排。"""
    configure_matplotlib_fonts()

    active_codes = sorted(set(active_codes))
    free_siso = sorted(set(free_siso))
    active_code_set = set(active_codes)
    matched_codes = set(final_match.keys())
    matched_siso = set(final_match.values())

    base_edge_color = "#b8b8b8"
    stage1_edge_color = "#7aa6c2"
    final_edge_color = "#2f4f4f"

    code_inactive_face = "white"
    code_active_wait_face = "#f2c879"
    code_active_match_face = "#74c69d"
    siso_idle_face = "white"
    siso_busy_face = "#6ea8d9"

    x_code = 0.10
    x_siso = 0.90
    rect_w = 0.14
    rect_pad = 0.02

    fig, ax = plt.subplots(figsize=(12, 9))
    fig.subplots_adjust(left=0.06, right=0.94, top=0.88, bottom=0.27)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    code_xy = [(x_code, y_pos(i, N_code)) for i in range(N_code)]
    siso_xy = [(x_siso, y_pos(j, N_siso)) for j in range(N_siso)]

    code_per_g = N_code // G
    siso_per_g = N_siso // G

    for g in range(G):
        c_start = g * code_per_g
        y0, h = group_rect_params(c_start, code_per_g, N_code)
        ax.add_patch(Rectangle((x_code - rect_w - rect_pad, y0), rect_w, h, fill=False, linewidth=2))

        s_start = g * siso_per_g
        y0s, hs = group_rect_params(s_start, siso_per_g, N_siso)
        ax.add_patch(Rectangle((x_siso + rect_pad, y0s), rect_w, hs, fill=False, linewidth=2))

        ax.text(x_code - rect_w - rect_pad - 0.03, y0 + h / 2, f"G{g}", va="center", ha="right", fontsize=11)
        ax.text(x_siso + rect_pad + rect_w + 0.03, y0s + hs / 2, f"G{g}", va="center", ha="left", fontsize=11)

    for code_idx, siso_idx in edges:
        x1, y1 = code_xy[code_idx]
        x2, y2 = siso_xy[siso_idx]
        ax.plot([x1, x2], [y1, y2], linewidth=0.8, alpha=0.25, color=base_edge_color, zorder=1)

    for code_idx, siso_idx in sorted(stage1_match.items()):
        x1, y1 = code_xy[code_idx]
        x2, y2 = siso_xy[siso_idx]
        ax.plot(
            [x1, x2],
            [y1, y2],
            linewidth=2.0,
            alpha=0.85,
            color=stage1_edge_color,
            linestyle="--",
            zorder=3,
        )

    for code_idx, siso_idx in sorted(final_match.items()):
        x1, y1 = code_xy[code_idx]
        x2, y2 = siso_xy[siso_idx]
        ax.plot([x1, x2], [y1, y2], linewidth=2.8, alpha=0.95, color=final_edge_color, zorder=4)

    for code_idx, (x, y) in enumerate(code_xy):
        if code_idx in matched_codes:
            face_color = code_active_match_face
            edge_width = 1.8
        elif code_idx in active_code_set:
            face_color = code_active_wait_face
            edge_width = 1.6
        else:
            face_color = code_inactive_face
            edge_width = 1.2

        ax.scatter(x, y, s=60, facecolors=face_color, edgecolors="black", linewidths=edge_width, zorder=5)
        font_weight = "bold" if code_idx in active_code_set else "normal"
        ax.text(x - 0.02, y, f"C{code_idx + 1}", va="center", ha="right", fontsize=9, fontweight=font_weight)

    free_siso_set = set(free_siso)
    for siso_idx, (x, y) in enumerate(siso_xy):
        if siso_idx in matched_siso:
            face_color = siso_busy_face
            edge_width = 1.8
        elif siso_idx in free_siso_set:
            face_color = siso_idle_face
            edge_width = 1.4
        else:
            face_color = "#dddddd"
            edge_width = 1.2

        ax.scatter(x, y, s=70, facecolors=face_color, edgecolors="black", linewidths=edge_width, zorder=5)
        ax.text(x + 0.02, y, f"S{siso_idx + 1}", va="center", ha="left", fontsize=9)

    max_possible = min(len(active_codes), len(free_siso))
    ax.text(
        0.5,
        1.03,
        (
            "Scheme B: Two-Stage Reconfigurable Scheduling"
            f" | Stage 1={len(stage1_match)}"
            f" | Final={len(final_match)}/{max_possible}"
        ),
        ha="center",
        va="bottom",
        transform=ax.transAxes,
        fontsize=13,
    )

    legend_items = [
        Line2D([0], [0], color=base_edge_color, lw=1.2, label="Allowed edges (topology background)"),
        Line2D([0], [0], color=stage1_edge_color, lw=2.0, linestyle="--", label="Stage 1: local initial matching"),
        Line2D([0], [0], color=final_edge_color, lw=2.8, label="Stage 2: final matching"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=code_active_match_face, markersize=8, lw=0, label="Active and finally matched code"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=code_active_wait_face, markersize=8, lw=0, label="Active but finally unmatched code"),
    ]
    ax.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False, fontsize=10)

    active_lines = chunk_labels("C", active_codes, chunk_size=8)
    waiting_lines = chunk_labels("C", waiting_codes, chunk_size=8)
    stage1_lines = format_matching_lines(stage1_match, chunk_size=4)
    final_lines = format_matching_lines(final_match, chunk_size=4)

    info_lines = [
        "Scheduling input:",
        f"  Active codes: {active_lines[0]}",
    ]
    info_lines.extend([f"             {line}" for line in active_lines[1:]])
    info_lines.append(f"  Stage 1 match count: {len(stage1_match)}")
    info_lines.append(f"  Stage 1 result: {stage1_lines[0]}")
    info_lines.extend([f"            {line}" for line in stage1_lines[1:]])
    info_lines.append(f"  Final match count: {len(final_match)}")
    info_lines.append(f"  Final result: {final_lines[0]}")
    info_lines.extend([f"            {line}" for line in final_lines[1:]])
    info_lines.append(f"  Waiting queue: {chunk_labels('C', waiting_codes, chunk_size=8)[0]}")
    waiting_lines = chunk_labels("C", waiting_codes, chunk_size=8)
    info_lines[-1] = f"  Waiting queue: {waiting_lines[0]}"
    info_lines.extend([f"           {line}" for line in waiting_lines[1:]])

    fig.text(0.07, 0.04, "\n".join(info_lines), ha="left", va="bottom", fontsize=10)

    if save:
        plt.savefig(save, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
