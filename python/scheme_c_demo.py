#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Visualize scheme C staged scheduling.

Outputs two figures by default:
- phase 1: normal codes scheduled inside each group
- phase 2: final schedule after tail-code local/bypass attempts
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Set, Tuple

from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt

from group_plotter import (
    chunk_labels,
    configure_matplotlib_fonts,
    format_matching_lines,
    group_rect_params,
    y_pos,
)
from group_scheduler import (
    DEFAULT_BYPASS_SCHEME_ID,
    get_bypass_edges,
    run_scheme_c_with_trace,
)


PLOT_BYPASS_SCHEME_LABELS = {
    "bypass_scheme_1": "Bypass Scheme 1",
    "bypass_scheme_2": "Bypass Scheme 2",
}


@dataclass
class Scenario:
    name: str
    description: str
    active_codes_1based: List[int]


def _tail_codes_for_group(n_code: int, g: int, group_g: int) -> List[int]:
    code_per_g = n_code // group_g
    code_begin = g * code_per_g
    return [code_begin + code_per_g - 2, code_begin + code_per_g - 1]


def _normal_codes_for_group(n_code: int, g: int, group_g: int) -> List[int]:
    code_per_g = n_code // group_g
    code_begin = g * code_per_g
    return list(range(code_begin, code_begin + code_per_g - 2))


def _bypass_edge_set(
    n_code: int,
    n_siso: int,
    extra_bypass_edges: Sequence[Tuple[int, int]],
) -> Set[Tuple[int, int]]:
    return {
        (code_idx, siso_idx)
        for code_idx, siso_idx in extra_bypass_edges
        if 0 <= code_idx < n_code and 0 <= siso_idx < n_siso
    }


def _chunk_labels_en(prefix: str, indices: Sequence[int], chunk_size: int = 8) -> List[str]:
    labels = [f"{prefix}{idx + 1}" for idx in indices]
    if not labels:
        return ["None"]
    return [", ".join(labels[i:i + chunk_size]) for i in range(0, len(labels), chunk_size)]


def _draw_scheme_c_stage(
    *,
    title: str,
    subtitle: str,
    n_code: int,
    n_siso: int,
    g: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
    background_edges: Sequence[Tuple[int, int]],
    match: Dict[int, int],
    waiting_codes: Sequence[int],
    save: str | None,
    show: bool,
    emphasize_bypass: bool,
    extra_bypass_edges: Sequence[Tuple[int, int]],
) -> None:
    configure_matplotlib_fonts()

    active_codes = sorted(set(active_codes))
    free_siso = sorted(set(free_siso))
    active_set = set(active_codes)
    matched_codes = set(match.keys())
    matched_siso = set(match.values())
    waiting_set = set(waiting_codes)
    bypass_edges = _bypass_edge_set(n_code, n_siso, extra_bypass_edges)

    x_code = 0.10
    x_siso = 0.90
    rect_w = 0.14
    rect_pad = 0.02

    base_edge_color = "#c9c9c9"
    bypass_edge_color = "#d28b36"
    match_edge_color = "#2f4f4f"

    code_inactive_face = "white"
    code_normal_face = "#8fd3c1"
    code_tail_face = "#f7c87c"
    code_wait_face = "#ef9a9a"
    siso_idle_face = "white"
    siso_busy_face = "#6ea8d9"

    fig, ax = plt.subplots(figsize=(12, 9))
    fig.subplots_adjust(left=0.06, right=0.94, top=0.88, bottom=0.26)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    code_xy = [(x_code, y_pos(i, n_code)) for i in range(n_code)]
    siso_xy = [(x_siso, y_pos(j, n_siso)) for j in range(n_siso)]

    code_per_g = n_code // g
    siso_per_g = n_siso // g

    tail_code_set: Set[int] = set()
    for group_idx in range(g):
        tail_code_set.update(_tail_codes_for_group(n_code, group_idx, g))

    for group_idx in range(g):
        c_start = group_idx * code_per_g
        y0, h = group_rect_params(c_start, code_per_g, n_code)
        ax.add_patch(Rectangle((x_code - rect_w - rect_pad, y0), rect_w, h, fill=False, linewidth=2))

        s_start = group_idx * siso_per_g
        y0s, hs = group_rect_params(s_start, siso_per_g, n_siso)
        ax.add_patch(Rectangle((x_siso + rect_pad, y0s), rect_w, hs, fill=False, linewidth=2))

        ax.text(x_code - rect_w - rect_pad - 0.03, y0 + h / 2, f"G{group_idx}", va="center", ha="right", fontsize=11)
        ax.text(x_siso + rect_pad + rect_w + 0.03, y0s + hs / 2, f"G{group_idx}", va="center", ha="left", fontsize=11)

    for code_idx, siso_idx in background_edges:
        x1, y1 = code_xy[code_idx]
        x2, y2 = siso_xy[siso_idx]
        is_bypass = (code_idx, siso_idx) in bypass_edges
        edge_color = bypass_edge_color if (emphasize_bypass and is_bypass) else base_edge_color
        edge_alpha = 0.35 if (emphasize_bypass and is_bypass) else 0.22
        edge_width = 1.0 if (emphasize_bypass and is_bypass) else 0.8
        ax.plot([x1, x2], [y1, y2], linewidth=edge_width, alpha=edge_alpha, color=edge_color, zorder=1)

    for code_idx, siso_idx in sorted(match.items()):
        x1, y1 = code_xy[code_idx]
        x2, y2 = siso_xy[siso_idx]
        ax.plot([x1, x2], [y1, y2], linewidth=2.8, alpha=0.95, color=match_edge_color, zorder=4)

    for code_idx, (x, y) in enumerate(code_xy):
        if code_idx in matched_codes:
            face_color = code_tail_face if code_idx in tail_code_set else code_normal_face
            edge_width = 1.8
        elif code_idx in waiting_set:
            face_color = code_wait_face
            edge_width = 1.8
        elif code_idx in active_set:
            face_color = code_tail_face if code_idx in tail_code_set else code_normal_face
            edge_width = 1.5
        else:
            face_color = code_inactive_face
            edge_width = 1.2

        ax.scatter(x, y, s=60, facecolors=face_color, edgecolors="black", linewidths=edge_width, zorder=5)
        font_weight = "bold" if code_idx in active_set else "normal"
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

    ax.text(0.5, 1.03, title, ha="center", va="bottom", transform=ax.transAxes, fontsize=13)
    ax.text(0.5, 0.995, subtitle, ha="center", va="top", transform=ax.transAxes, fontsize=10)

    legend_items = [
        Line2D([0], [0], color=base_edge_color, lw=1.2, label="Local allowed edges"),
        Line2D([0], [0], color=bypass_edge_color, lw=1.2, label="Bypass edges"),
        Line2D([0], [0], color=match_edge_color, lw=2.8, label="Scheduled edges in this stage"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=code_normal_face, markersize=8, lw=0, label="Normal code"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=code_tail_face, markersize=8, lw=0, label="Tail code"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor=code_wait_face, markersize=8, lw=0, label="Active but unscheduled"),
    ]
    ax.legend(handles=legend_items, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=3, frameon=False, fontsize=10)

    active_lines = _chunk_labels_en("C", active_codes, chunk_size=8)
    waiting_lines = _chunk_labels_en("C", waiting_codes, chunk_size=8)
    match_lines = format_matching_lines(match, chunk_size=4)

    info_lines = [
        "Scheduling input:",
        f"  Active codes: {active_lines[0]}",
    ]
    info_lines.extend([f"                {line}" for line in active_lines[1:]])
    info_lines.append(f"  Scheduled count: {len(match)}")
    info_lines.append(f"  Match result: {match_lines[0]}")
    info_lines.extend([f"               {line}" for line in match_lines[1:]])
    info_lines.append(f"  Unscheduled: {waiting_lines[0]}")
    info_lines.extend([f"               {line}" for line in waiting_lines[1:]])

    fig.text(0.07, 0.04, "\n".join(info_lines), ha="left", va="bottom", fontsize=10)

    if save:
        plt.savefig(save, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def main() -> None:
    n_code = 32
    n_siso = 16
    g = 4
    free_siso_1based = list(range(1, n_siso + 1))
    show_figure = False
    scheme_ids = [DEFAULT_BYPASS_SCHEME_ID, "bypass_scheme_2"]

    free_siso = [idx - 1 for idx in free_siso_1based]
    scenarios = [
        Scenario(
            name="low_load",
            description="Fewer active codes; most groups are not saturated.",
            active_codes_1based=[1, 3, 7, 11, 15, 19, 23, 31],
        ),
        Scenario(
            name="high_load",
            description="More active codes; local SISOs are likely to saturate before tail-code scheduling.",
            active_codes_1based=[1, 2, 3, 4, 5, 7, 8, 11, 12, 13, 15, 16, 19, 23, 24, 31, 32],
        ),
    ]

    for scheme_id in scheme_ids:
        scheme_label = PLOT_BYPASS_SCHEME_LABELS.get(scheme_id, scheme_id)
        extra_bypass_edges = get_bypass_edges(scheme_id)

        for scenario in scenarios:
            active_codes = [idx - 1 for idx in scenario.active_codes_1based]
            stage1_match, final_match, waiting_codes, local_edges, all_edges = run_scheme_c_with_trace(
                N_code=n_code,
                N_siso=n_siso,
                G=g,
                active_codes=active_codes,
                free_siso=free_siso,
                extra_bypass_edges=extra_bypass_edges,
            )
            stage1_waiting = [code_idx for code_idx in active_codes if code_idx not in stage1_match]

            phase1_save_path = f"schemeC_{scenario.name}_phase1_4group_{scheme_id}.png"
            phase2_save_path = f"schemeC_{scenario.name}_phase2_4group_{scheme_id}.png"

            _draw_scheme_c_stage(
                title=f"Scheme C Phase 1 ({scenario.name}, {scheme_label}) | scheduled={len(stage1_match)}",
                subtitle=scenario.description + " Only local SISOs are used in phase 1.",
                n_code=n_code,
                n_siso=n_siso,
                g=g,
                active_codes=active_codes,
                free_siso=free_siso,
                background_edges=local_edges,
                match=stage1_match,
                waiting_codes=stage1_waiting,
                save=phase1_save_path,
                show=show_figure,
                emphasize_bypass=False,
                extra_bypass_edges=extra_bypass_edges,
            )

            _draw_scheme_c_stage(
                title=f"Scheme C Phase 2 ({scenario.name}, {scheme_label}) | scheduled={len(final_match)}",
                subtitle=scenario.description + " Tail codes can use local free SISOs and bypass edges in phase 2.",
                n_code=n_code,
                n_siso=n_siso,
                g=g,
                active_codes=active_codes,
                free_siso=free_siso,
                background_edges=all_edges,
                match=final_match,
                waiting_codes=waiting_codes,
                save=phase2_save_path,
                show=show_figure,
                emphasize_bypass=True,
                extra_bypass_edges=extra_bypass_edges,
            )

            print(f"Scheme C demo completed for scenario: {scenario.name}")
            print(f"Bypass scheme: {scheme_label}")
            print(f"Description: {scenario.description}")
            print("Active Codes:", [f"C{i + 1}" for i in active_codes])
            print("Phase 1:", [f"C{code_idx + 1}->S{siso_idx + 1}" for code_idx, siso_idx in sorted(stage1_match.items())])
            print("Final   :", [f"C{code_idx + 1}->S{siso_idx + 1}" for code_idx, siso_idx in sorted(final_match.items())])
            print("Waiting :", [f"C{i + 1}" for i in waiting_codes])
            print(f"Saved figure to: {phase1_save_path}")
            print(f"Saved figure to: {phase2_save_path}")
            print()


if __name__ == "__main__":
    main()
