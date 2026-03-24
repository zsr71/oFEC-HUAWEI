#pragma once

namespace newcode {
namespace detail {

// ----- find L least reliable core positions (indices 0..254) by |LLR| -----
template<typename LLR>
static void find_least_reliable(const LLR* LLR_list, int L,
                                std::vector<int>& pos, std::vector<float>& absval)
{
    // 输入:
    // - LLR_list: 长度至少为 BCH_N_CORE 的输入 LLR 序列。
    // - L: 需要保留的最不可靠位数量。
    // - pos/absval: 输出容器引用。
    // 输出:
    // - pos: 绝对值最小的若干位置。
    // - absval: 对应位置的 |LLR|。
    // 用途:
    // - Chase 首先需要挑出最不可靠位，再只对这些位做翻转测试。
    struct Node { float a; int i; };
    std::vector<Node> v; v.reserve(BCH_N_CORE);
    for (int i = 0; i < BCH_N_CORE; ++i)
        v.push_back({ std::fabs(llr_to_float(LLR_list[i])), i });

    const int take = std::min(L, (int)v.size());
    if (take < (int)v.size())
        std::nth_element(v.begin(), v.begin() + take, v.end(),
                         [](const Node& x, const Node& y){ return x.a < y.a; });
    std::sort(v.begin(), v.begin() + take, [](const Node& x, const Node& y){ return x.a < y.a; });

    pos.resize(take);
    absval.resize(take);
    for (int k = 0; k < take; ++k) { pos[k] = v[k].i; absval[k] = v[k].a; }
}

} // namespace detail
} // namespace newcode
