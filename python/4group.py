#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4group.py
运行入口：
- 定义实验参数
- 调用调度逻辑
- 调用画图模块生成方案 A / B 图
"""

from group_plotter import draw_scheme_a_schedule, draw_scheme_b_schedule
from group_scheduler import run_scheme_a, run_scheme_b


def main() -> None:
    # -----------------------------
    # 这里直接定义当前实验配置。
    # 如果要改测试场景，直接修改下面这些变量即可。
    # 编号统一使用 1-based，更贴近图上的 C1/S1 标注。
    # -----------------------------
    N_code = 32
    N_siso = 16
    G = 4
    active_codes_1based = [1, 2, 3, 4, 5,7, 8, 11, 12, 13, 15, 16, 19, 23, 24, 31, 32]
    free_siso_1based = list(range(1, N_siso + 1))
    show_figure = True

    active_codes = [idx - 1 for idx in active_codes_1based]
    free_siso = [idx - 1 for idx in free_siso_1based]

    scheme_a_save_path = "schemeA_4group.png"
    scheme_b_save_path = "schemeB_4group.png"

    match_a, waiting_codes_a, all_edges_a = run_scheme_a(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    draw_scheme_a_schedule(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
        edges=all_edges_a,
        match=match_a,
        waiting_codes=waiting_codes_a,
        save=scheme_a_save_path,
        show=show_figure,
    )
    print("方案A调度完成。")
    print("活跃 Code:", [f"C{i + 1}" for i in active_codes])
    print("匹配结果:", [f"C{code_idx + 1}->S{siso_idx + 1}" for code_idx, siso_idx in sorted(match_a.items())])
    print("等待队列:", [f"C{i + 1}" for i in waiting_codes_a])
    print(f"Saved figure to: {scheme_a_save_path}")

    stage1_match_b, final_match_b, waiting_codes_b, all_edges_b = run_scheme_b(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    draw_scheme_b_schedule(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
        edges=all_edges_b,
        stage1_match=stage1_match_b,
        final_match=final_match_b,
        waiting_codes=waiting_codes_b,
        save=scheme_b_save_path,
        show=show_figure,
    )
    print("方案B调度完成。")
    print("活跃 Code:", [f"C{i + 1}" for i in active_codes])
    print("阶段1结果:", [f"C{code_idx + 1}->S{siso_idx + 1}" for code_idx, siso_idx in sorted(stage1_match_b.items())])
    print("最终结果:", [f"C{code_idx + 1}->S{siso_idx + 1}" for code_idx, siso_idx in sorted(final_match_b.items())])
    print("等待队列:", [f"C{i + 1}" for i in waiting_codes_b])
    print(f"Saved figure to: {scheme_b_save_path}")


if __name__ == "__main__":
    main()
