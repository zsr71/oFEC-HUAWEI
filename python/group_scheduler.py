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


def _build_bypass_map(
    N_code: int,
    N_siso: int,
    extra_edges: Sequence[Tuple[int, int]],
) -> Dict[int, List[int]]:
    """把旁路线列表转成 code -> bypass siso 列表。"""
    bypass_map: Dict[int, List[int]] = {i: [] for i in range(N_code)}
    for code_idx, siso_idx in extra_edges:
        if 0 <= code_idx < N_code and 0 <= siso_idx < N_siso:
            bypass_map[code_idx].append(siso_idx)
    for code_idx in bypass_map:
        bypass_map[code_idx] = sorted(set(bypass_map[code_idx]))
    return bypass_map


def schedule_scheme_c_staged(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], List[int]]:
    """方案 C：返回最终结果的简化接口。"""
    stage1_match, final_match, waiting_codes = schedule_scheme_c_staged_with_trace(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    _ = stage1_match
    return final_match, waiting_codes


def schedule_scheme_c_staged_with_trace(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], Dict[int, int], List[int]]:
    """
    方案 C：顺序式两阶段调度。

    仅针对当前固定的 G=4 场景：
    - 每组 8 个 code
    - 每组 4 个 siso
    - 第一阶段先调度每组前 6 个普通 code
    - 第二阶段再调度每组后 2 个边界 code
      * 先尝试本组空闲 siso
      * 再尝试额外旁路线
    """
    code_per_g, siso_per_g = _validate_group_shape(N_code, N_siso, G)
    if code_per_g < 2:
        raise ValueError("每组至少要有 2 个 code 才能区分普通码字和边界码字")

    active_set = set(active_codes)
    free_siso_set = set(free_siso)
    code_to_siso: Dict[int, int] = {}
    used_siso: Set[int] = set()
    bypass_map = _build_bypass_map(N_code, N_siso, EXTRA_BYPASS_EDGES)

    def try_assign(code_idx: int, candidates: Sequence[int]) -> bool:
        for siso_idx in candidates:
            if siso_idx not in free_siso_set or siso_idx in used_siso:
                continue
            code_to_siso[code_idx] = siso_idx
            used_siso.add(siso_idx)
            return True
        return False

    normal_count = code_per_g - 2
    stage1_match: Dict[int, int] = {}

    # Phase 1: 每组前 6 个普通码字先调度，只用本组 SISO
    for g in range(G):
        code_begin = g * code_per_g
        local_sisos = list(range(g * siso_per_g, (g + 1) * siso_per_g))
        normal_codes = list(range(code_begin, code_begin + normal_count))

        for code_idx in normal_codes:
            if code_idx not in active_set:
                continue
            assigned = try_assign(code_idx, local_sisos)
            if assigned:
                stage1_match[code_idx] = code_to_siso[code_idx]
            if all(siso in used_siso for siso in local_sisos):
                break

    # Phase 2: 每组最后 2 个边界码字后调度，先本组后旁路
    for g in range(G):
        code_begin = g * code_per_g
        local_sisos = list(range(g * siso_per_g, (g + 1) * siso_per_g))
        tail_codes = list(range(code_begin + normal_count, code_begin + code_per_g))

        for code_idx in tail_codes:
            if code_idx not in active_set:
                continue
            if code_idx in code_to_siso:
                continue
            if try_assign(code_idx, local_sisos):
                continue
            try_assign(code_idx, bypass_map.get(code_idx, []))

    waiting_codes = [code_idx for code_idx in sorted(active_set) if code_idx not in code_to_siso]
    return stage1_match, code_to_siso, waiting_codes


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


def run_scheme_c(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], List[int], List[Tuple[int, int]]]:
    """方案 C：顺序式两阶段调度。"""
    active_codes = sorted(set(active_codes))
    free_siso = sorted(set(free_siso))
    all_edges = build_allowed_edges(N_code, N_siso, G)
    final_match, waiting_codes = schedule_scheme_c_staged(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    return final_match, waiting_codes, all_edges


def run_scheme_c_with_trace(
    N_code: int,
    N_siso: int,
    G: int,
    active_codes: Sequence[int],
    free_siso: Sequence[int],
) -> Tuple[Dict[int, int], Dict[int, int], List[int], List[Tuple[int, int]], List[Tuple[int, int]]]:
    """方案 C：返回阶段1结果、最终结果和可视化所需边集。"""
    active_codes = sorted(set(active_codes))
    free_siso = sorted(set(free_siso))
    local_edges = build_local_edges(N_code, N_siso, G)
    all_edges = build_allowed_edges(N_code, N_siso, G)
    stage1_match, final_match, waiting_codes = schedule_scheme_c_staged_with_trace(
        N_code=N_code,
        N_siso=N_siso,
        G=G,
        active_codes=active_codes,
        free_siso=free_siso,
    )
    return stage1_match, final_match, waiting_codes, local_edges, all_edges
