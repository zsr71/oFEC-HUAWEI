#pragma once

namespace newcode {
namespace detail {

// ----- find L least reliable core positions (indices 0..254) by |LLR| -----
template<typename LLR>
static void find_least_reliable(const LLR* LLR_list, int L,
                                std::vector<int>& pos, std::vector<float>& absval)
{
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

