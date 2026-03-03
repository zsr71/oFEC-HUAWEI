import matplotlib.pyplot as plt

NUM_RED = 32
NUM_BLUE = 16

def build_edges(N, num_red=32, num_blue=16):
    """
    按如下构造生成边：
    - 将 32 个红点分成 16 组，每组 2 个
    - 蓝点 B_i 连接到连续 N 组红点
    - 因此每个蓝点连接 2N 个红点，每个红点连接 N 个蓝点
    """
    if num_red != 2 * num_blue:
        raise ValueError("当前这套构造默认要求 num_red = 2 * num_blue，例如 32 和 16。")
    if not (0 <= N <= num_blue):
        raise ValueError(f"N 必须满足 0 <= N <= {num_blue}")

    edges = []

    # 蓝点编号用 0~15，对应 B1~B16
    # 红点分组编号用 0~15，对应:
    # group 0 -> R1,R2
    # group 1 -> R3,R4
    # ...
    # group 15 -> R31,R32
    for b in range(num_blue):
        for k in range(N):
            g = (b + k) % num_blue  # 连续 N 组，按模 num_blue 循环
            r1 = 2 * g + 1
            r2 = 2 * g + 2
            edges.append((r1, b + 1))
            edges.append((r2, b + 1))

    return edges


def check_degrees(edges, num_red=32, num_blue=16):
    """
    检查每个红点和蓝点的度数
    """
    red_deg = {i: 0 for i in range(1, num_red + 1)}
    blue_deg = {j: 0 for j in range(1, num_blue + 1)}

    for r, b in edges:
        red_deg[r] += 1
        blue_deg[b] += 1

    return red_deg, blue_deg


def plot_bipartite(N, ax=None, num_red=32, num_blue=16):
    """
    画出二分图连接方式
    """
    edges = build_edges(N, num_red, num_blue)
    red_deg, blue_deg = check_degrees(edges, num_red, num_blue)

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 12))

    # 左右两列点的位置
    x_red = 0
    x_blue = 1

    # 为了让 R1 在上方，这里把 y 轴倒过来
    red_y = {i: num_red - i + 1 for i in range(1, num_red + 1)}
    blue_y = {j: (num_blue - j) * (num_red - 1) / (num_blue - 1) + 1 for j in range(1, num_blue + 1)}

    # 画边
    for r, b in edges:
        ax.plot(
            [x_red, x_blue],
            [red_y[r], blue_y[b]],
            linewidth=0.8,
            alpha=0.45
        )

    # 画红点
    ax.scatter(
        [x_red] * num_red,
        [red_y[i] for i in range(1, num_red + 1)],
        s=60,
        c='red',
        label='Red nodes'
    )

    # 画蓝点
    ax.scatter(
        [x_blue] * num_blue,
        [blue_y[j] for j in range(1, num_blue + 1)],
        s=80,
        c='blue',
        label='Blue nodes'
    )

    # 标注红点编号
    for i in range(1, num_red + 1):
        ax.text(x_red - 0.03, red_y[i], f"R{i}", ha='right', va='center', fontsize=8)

    # 标注蓝点编号
    for j in range(1, num_blue + 1):
        ax.text(x_blue + 0.03, blue_y[j], f"B{j}", ha='left', va='center', fontsize=8)

    # show degree information in title
    red_deg_set = sorted(set(red_deg.values()))
    blue_deg_set = sorted(set(blue_deg.values()))
    ax.set_title(
        f"N = {N} | red degrees = {red_deg_set} | blue degrees = {blue_deg_set}",
        fontsize=12
    )

    ax.set_xlim(-0.2, 1.2)
    ax.set_ylim(0, num_red + 1)
    ax.axis('off')


def plot_multiple_N(N_list, num_red=32, num_blue=16):
    """
    一次显示多个不同 N 的连接图
    """
    cols = 2
    rows = (len(N_list) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(14, 6 * rows))
    axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]

    for ax, N in zip(axes, N_list):
        plot_bipartite(N, ax=ax, num_red=num_red, num_blue=num_blue)

    # 多余子图关掉
    for i in range(len(N_list), len(axes)):
        axes[i].axis('off')

    plt.tight_layout()
    plt.show()


def print_connection_table(N, num_red=32, num_blue=16):
    """
    打印每个蓝点连接哪些红点
    """
    edges = build_edges(N, num_red, num_blue)
    blue_to_red = {b: [] for b in range(1, num_blue + 1)}

    for r, b in edges:
        blue_to_red[b].append(r)

    print(f"\n===== N = {N} 的连接方案 =====")
    for b in range(1, num_blue + 1):
        red_list = sorted(blue_to_red[b])
        print(f"B{b:>2} -> {red_list}")


if __name__ == "__main__":
    # 1) 画单个 N
    plt.figure(figsize=(10, 12))
    plot_bipartite(N=3)
    plt.show()

    # 2) 一次画多个不同 N
    plot_multiple_N([1, 2, 3, 4, 6, 8])

    # 3) 打印某个 N 的具体编号方案
    print_connection_table(3)