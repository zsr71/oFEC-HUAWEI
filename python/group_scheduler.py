#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_scheduler.py
4 组基本隔离 + 少量跨组旁路的拓扑与调度逻辑。
这里不包含任何 Matplotlib 依赖。
"""

from typing import Dict, List, Sequence, Set, Tuple


# 额外开放的跨组旁路线，采用 0-based 索引。
EXTRA_BYPASS_EDGES = [
    (6, 11),   # C7  -> S12
    (6, 15),   # C7  -> S16
    (7, 11),   # C8  -> S12
    (7, 15),   # C8  -> S16
    (14, 11),  # C15 -> S12
    (14, 15),  # C15 -> S16
    (15, 11),  # C16 -> S12
    (15, 15),  # C16 -> S16
    (22, 3),   # C23 -> S4
    (22, 7),   # C23 -> S8
    (23, 3),   # C24 -> S4
    (23, 7),   # C24 -> S8
    (30, 3),   # C31 -> S4
    (30, 7),   # C31 -> S8
    (31, 3),   # C32 -> S4
    (31, 7),   # C32 -> S8
]


def _validate_group_shape(N_code: int, N_siso: int, G: int) -> Tuple[int, int]:
    """检查组划分是否合法，并返回每组的 code / siso 数量。"""
    if N_code % G != 0 or N_siso % G != 0:
        raise ValueError(
            f"N_code 和 N_siso 必须能被 G 整除，当前 N_code={N_code}, N_siso={N_siso}, G={G}"
        )
    return N_code // G, N_siso // G


def build_allowed_edges(N_code: int, N_siso: int, G: int) -> List[Tuple[int, int]]:
    """构造当前拓扑允许的全部边：组内边 + 额外旁路线。"""
    code_per_g, siso_per_g = _validate_group_shape(N_code, N_siso, G)
    edges: Set[Tuple[int, int]] = set()

    for i in range(N_code):
        g = i // code_per_g
        s_start = g * siso_per_g
        s_end = (g + 1) * siso_per_g
        for j in range(s_start, s_end):
            edges.add((i, j))

    for code_idx, siso_idx in EXTRA_BYPASS_EDGES:
        if 0 <= code_idx < N_code and 0 <= siso_idx < N_siso:
            edges.add((code_idx, siso_idx))

    return sorted(edges)


def build_local_edges(N_code: int, N_siso: int, G: int) -> List[Tuple[int, int]]:
    """只构造组内边，不包含任何跨组旁路线。"""
    code_per_g, siso_per_g = _validate_group_shape(N_code, N_siso, G)
    edges: Set[Tuple[int, int]] = set()

    for i in range(N_code):
        g = i // code_per_g
        s_start = g * siso_per_g
        s_end = (g + 1) * siso_per_g
        for j in range(s_start, s_end):
            edges.add((i, j))

    return sorted(edges)


def build_adjacency(edges: Sequence[Tuple[int, int]], N_code: int) -> Dict[int, List[int]]:
    """把边列表转成邻接表，便于后续做最大匹配。"""
    adjacency = {i: [] for i in range(N_code)}
    for code_idx, siso_idx in edges:
        adjacency[code_idx].append(siso_idx)
    return adjacency


def try_augment(
    code_idx: int,
    adjacency: Dict[int, List[int]],
    free_siso: Set[int],
    match_siso_to_code: Dict[int, int],
    visited_siso: Set[int],
) -> bool:
    """
    DFS 增广：
    尝试为当前 code 找到一个可用的 SISO；
    如果目标 SISO 已被占用，就递归尝试把原 holder 挪走。
    """
    for siso_idx in adjacency[code_idx]:
        if siso_idx not in free_siso or siso_idx in visited_siso:
            continue

        visited_siso.add(siso_idx)
        holder = match_siso_to_code.get(siso_idx)

        if holder is None or try_augment(holder, adjacency, free_siso, match_siso_to_code, visited_siso):
            match_siso_to_code[siso_idx] = code_idx
            return True

    return False


def maximum_bipartite_matching(
    active_codes: Sequence[int],
    adjacency: Dict[int, List[int]],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], List[int]]:
    """
    对当前活跃 code 与空闲 SISO 做一次最大匹配。

    返回：
    - code_to_siso: 已匹配成功的 code -> siso
    - waiting_codes: 本轮仍未匹配到的活跃 code
    """
    free_siso_set = set(free_siso)
    match_siso_to_code: Dict[int, int] = {}

    for code_idx in active_codes:
        visited_siso: Set[int] = set()
        try_augment(code_idx, adjacency, free_siso_set, match_siso_to_code, visited_siso)

    code_to_siso = {code_idx: siso_idx for siso_idx, code_idx in match_siso_to_code.items()}
    waiting_codes = [code_idx for code_idx in active_codes if code_idx not in code_to_siso]
    return code_to_siso, waiting_codes


def augment_from_seed_matching(
    active_codes: Sequence[int],
    adjacency: Dict[int, List[int]],
    free_siso: Sequence[int],
    seed_code_to_siso: Dict[int, int],
) -> Tuple[Dict[int, int], List[int]]:
    """
    在已有匹配基础上继续增广。

    这个函数对应方案 B 的“阶段 2 可重排版”：
    - 先拿阶段 1 的组内匹配作为初始解
    - 再开放全边集（组内边 + 旁路线）
    - 允许沿增广路重排已有匹配，以便继续提高匹配数
    """
    free_siso_set = set(free_siso)
    match_siso_to_code = {siso_idx: code_idx for code_idx, siso_idx in seed_code_to_siso.items()}

    for code_idx in active_codes:
        if code_idx in seed_code_to_siso:
            continue
        visited_siso: Set[int] = set()
        try_augment(code_idx, adjacency, free_siso_set, match_siso_to_code, visited_siso)

    code_to_siso = {code_idx: siso_idx for siso_idx, code_idx in match_siso_to_code.items()}
    waiting_codes = [code_idx for code_idx in active_codes if code_idx not in code_to_siso]
    return code_to_siso, waiting_codes


def schedule_scheme_b_reconfig(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], Dict[int, int], List[int]]:
    """
    方案 B：两阶段 + 可重排。

    返回：
    - stage1_match: 仅用组内边得到的初始匹配
    - final_match: 开放旁路后，在初始匹配基础上继续增广得到的最终匹配
    - waiting_codes: 最终仍未匹配的活跃 code
    """
    local_edges = build_local_edges(N_code, N_siso, G)
    all_edges = build_allowed_edges(N_code, N_siso, G)

    local_adjacency = build_adjacency(local_edges, N_code)
    all_adjacency = build_adjacency(all_edges, N_code)

    stage1_match, _ = maximum_bipartite_matching(active_codes, local_adjacency, free_siso)
    final_match, waiting_codes = augment_from_seed_matching(active_codes, all_adjacency, free_siso, stage1_match)
    return stage1_match, final_match, waiting_codes


def run_scheme_a(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], List[int], List[Tuple[int, int]]]:
    """方案 A：全局最大匹配。"""
    active_codes = sorted(set(active_codes))
    free_siso = sorted(set(free_siso))
    all_edges = build_allowed_edges(N_code, N_siso, G)
    adjacency = build_adjacency(all_edges, N_code)
    match, waiting_codes = maximum_bipartite_matching(active_codes, adjacency, free_siso)
    return match, waiting_codes, all_edges


def run_scheme_b(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], Dict[int, int], List[int], List[Tuple[int, int]]]:
    """方案 B：两阶段 + 可重排。"""
    active_codes = sorted(set(active_codes))
    free_siso = sorted(set(free_siso))
    all_edges = build_allowed_edges(N_code, N_siso, G)
    stage1_match, final_match, waiting_codes = schedule_scheme_b_reconfig(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    return stage1_match, final_match, waiting_codes, all_edges
