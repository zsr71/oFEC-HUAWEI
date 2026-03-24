#pragma once

#include "newcode/common/qfloat/qfloat.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace qfloat {

template<int NBITS, typename Store>
constexpr int qfloat<NBITS, Store>::Q()
{
    // 输入: 无。
    // 输出: 正饱和码值上界 Q，NBITS=5 时返回 15。
    // 用途: 统一定义量化码值范围 [-Q, +Q]。
    return (1 << (NBITS - 1)) - 1;
}

template<int NBITS, typename Store>
constexpr int qfloat<NBITS, Store>::LO()
{
    // 输入: 无。
    // 输出: 负饱和码值下界，等于 -Q()。
    // 用途: 用于饱和裁剪和边界判断。
    return -Q();
}

template<int NBITS, typename Store>
constexpr int qfloat<NBITS, Store>::HI()
{
    // 输入: 无。
    // 输出: 正饱和码值上界，等于 +Q()。
    // 用途: 用于饱和裁剪和边界判断。
    return +Q();
}

template<int NBITS, typename Store>
float qfloat<NBITS, Store>::current_clip()
{
    // 输入: 无。
    // 输出: 当前线程本地的 clip 浮点值。
    // 用途: 控制 float <-> qfloat 转换时的量化满量程。
    return current_clip_ref();
}

template<int NBITS, typename Store>
void qfloat<NBITS, Store>::set_clip(float clip)
{
    // 输入: clip，量化时使用的满量程幅度。
    // 输出: 无。
    // 用途: 在一段处理流程内统一量化/反量化尺度。
    current_clip_ref() = clip;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>::qfloat() : code_(0)
{
    // 输入: 无。
    // 输出: 一个内部 code_=0 的对象。
    // 用途: 表示 0 LLR 或作为容器默认值。
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>::qfloat(int code, float clip) : code_(sat(code))
{
    // 输入: code 为目标码值，clip 参数在这里不参与计算。
    // 输出: 一个内部码值为 sat(code) 的对象。
    // 用途: 当上层已经持有码值时，直接恢复 qfloat 表示。
    (void)clip;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>::qfloat(float x, float clip) : code_(0)
{
    // 输入: x 为真实浮点值，clip 为量化满量程。
    // 输出: 一个内部码值为 quantize(x, clip) 的对象。
    // 用途: 把真实 LLR 或其他浮点量映射到线性量化码值。
    code_ = quantize(x, clip);
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> qfloat<NBITS, Store>::from_float(float x, float clip)
{
    // 输入: x 为待量化浮点值，clip 为量化满量程。
    // 输出: 对应的 qfloat 对象。
    // 用途: 比显式调用构造函数更清晰地表达“从 float 量化”。
    return qfloat<NBITS, Store>(x, clip);
}

template<int NBITS, typename Store>
float qfloat<NBITS, Store>::to_float() const
{
    // 输入: 当前对象内部的 code_，以及当前线程 clip。
    // 输出: 近似还原后的 float 值。
    // 用途: 在需要真实幅度时，把量化表示还原到浮点域。
    const float clip = current_clip();
    return static_cast<float>(code_) * (clip / static_cast<float>(Q()));
}

template<int NBITS, typename Store>
float qfloat<NBITS, Store>::clip()
{
    // 输入: 无。
    // 输出: 当前 clip 浮点值。
    // 用途: 方便按类型接口访问 clip 参数。
    return current_clip();
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> qfloat<NBITS, Store>::operator-() const
{
    // 输入: 当前对象的 code_。
    // 输出: 码值为 sat(-code_) 的新对象。
    // 用途: 保持量化域内的符号反转，并自动饱和。
    return qfloat<NBITS, Store>(sat(-code_));
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator+=(qfloat rhs)
{
    // 输入: rhs 为另一个 qfloat。
    // 输出: 返回 *this，内部 code_ 更新为饱和后的 code_ + rhs.code_。
    // 用途: 在码值域近似执行线性加法，例如 LLR 合并。
    code_ = sat(int(code_) + int(rhs.code_));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator-=(qfloat rhs)
{
    // 输入: rhs 为另一个 qfloat。
    // 输出: 返回 *this，内部 code_ 更新为饱和后的 code_ - rhs.code_。
    // 用途: 在码值域近似执行线性减法。
    code_ = sat(int(code_) - int(rhs.code_));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator*=(float k)
{
    // 输入: k 为缩放系数。
    // 输出: 返回 *this，内部 code_ 按 code_*k 四舍五入并饱和。
    // 用途: 对量化值做增益缩放，例如 alpha/beta 类参数调整。
    code_ = sat(static_cast<int>(std::lrint(static_cast<float>(code_) * k)));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>& qfloat<NBITS, Store>::operator/=(float k)
{
    // 输入: k 为除数；若 k 为 0，则保持原值不变。
    // 输出: 返回 *this，内部 code_ 按 code_/k 四舍五入并饱和。
    // 用途: 在码值域做比例缩放，同时避免除零异常。
    if (k != 0.0f)
        code_ = sat(static_cast<int>(std::lrint(static_cast<float>(code_) / k)));
    return *this;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator+(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 饱和后的和。
    // 用途: 提供与内置数值类型一致的 a + b 写法。
    a += b;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator-(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 饱和后的差。
    // 用途: 提供与内置数值类型一致的 a - b 写法。
    a -= b;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(qfloat<NBITS, Store> a, float k)
{
    // 输入: a 为量化值，k 为缩放系数。
    // 输出: 缩放并饱和后的 qfloat。
    // 用途: 支持 a * k 写法。
    a *= k;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(float k, qfloat<NBITS, Store> a)
{
    // 输入: k 为缩放系数，a 为量化值。
    // 输出: 缩放并饱和后的 qfloat。
    // 用途: 支持 k * a 写法。
    a *= k;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator/(qfloat<NBITS, Store> a, float k)
{
    // 输入: a 为量化值，k 为除数。
    // 输出: 缩放并饱和后的 qfloat。
    // 用途: 支持 a / k 写法。
    a /= k;
    return a;
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator*(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: a、b 为两个量化值。
    // 输出: 先反量化到 float 相乘，再重新量化后的 qfloat。
    // 用途: 对乘法保留更接近真实幅度的结果，而不是简单相乘码值。
    float y = a.to_float() * b.to_float();
    return qfloat<NBITS, Store>::from_float(y);
}

template<int NBITS, typename Store>
qfloat<NBITS, Store> operator/(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: a 为被除数，b 为除数。
    // 输出: 先反量化到 float 相除，再重新量化后的 qfloat；若分母为 0，直接返回 a。
    // 用途: 对除法保留真实数值语义，并避免除零。
    float den = b.to_float();
    if (den == 0.0f) return a;
    float y = a.to_float() / den;
    return qfloat<NBITS, Store>::from_float(y);
}

template<int NBITS, typename Store>
bool operator<(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 比较结果布尔值。
    // 用途: 在线性量化坐标系中做顺序关系判断。
    return a.code() < b.code();
}

template<int NBITS, typename Store>
bool operator>(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 比较结果布尔值。
    // 用途: 在线性量化坐标系中做顺序关系判断。
    return a.code() > b.code();
}

template<int NBITS, typename Store>
bool operator<=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 比较结果布尔值。
    // 用途: 在线性量化坐标系中做顺序关系判断。
    return a.code() <= b.code();
}

template<int NBITS, typename Store>
bool operator>=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 比较结果布尔值。
    // 用途: 在线性量化坐标系中做顺序关系判断。
    return a.code() >= b.code();
}

template<int NBITS, typename Store>
bool operator==(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 比较结果布尔值。
    // 用途: 判断两个量化对象的内部码值是否完全一致。
    return a.code() == b.code();
}

template<int NBITS, typename Store>
bool operator!=(qfloat<NBITS, Store> a, qfloat<NBITS, Store> b)
{
    // 输入: 两个 qfloat。
    // 输出: 比较结果布尔值。
    // 用途: 判断两个量化对象的内部码值是否不同。
    return a.code() != b.code();
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>::operator float() const
{
    // 输入: 当前 qfloat 对象。
    // 输出: 反量化后的 float。
    // 用途: 需要真实幅度时可直接 static_cast<float>(q)。
    return to_float();
}

template<int NBITS, typename Store>
qfloat<NBITS, Store>::operator int() const
{
    // 输入: 当前 qfloat 对象。
    // 输出: 内部量化码值。
    // 用途: 需要直接读取码值时可 static_cast<int>(q)。
    return static_cast<int>(code_);
}

template<int NBITS, typename Store>
int qfloat<NBITS, Store>::code() const
{
    // 输入: 当前 qfloat 对象。
    // 输出: int 类型的 code_。
    // 用途: 在码值域算法里直接使用量化索引或整数和。
    return static_cast<int>(code_);
}

template<int NBITS, typename Store>
void qfloat<NBITS, Store>::set_code(int c)
{
    // 输入: c 为目标码值。
    // 输出: 无，内部 code_ 更新为 sat(c)。
    // 用途: 当上层已经在码值域完成运算后，回填到 qfloat 对象。
    code_ = sat(c);
}

template<int NBITS, typename Store>
float& qfloat<NBITS, Store>::current_clip_ref()
{
    // 输入: 无。
    // 输出: thread_local float 引用。
    // 用途: 让不同线程可独立维护量化满量程。
    thread_local float clip = 0.0f;
    return clip;
}

template<int NBITS, typename Store>
int qfloat<NBITS, Store>::sat(int x)
{
    // 输入: x 为任意整数码值。
    // 输出: 被限制到 [LO(), HI()] 区间内的整数。
    // 用途: 保证内部 code_ 不越界。
    if (x > HI()) return HI();
    if (x < LO()) return LO();
    return x;
}

template<int NBITS, typename Store>
Store qfloat<NBITS, Store>::quantize(float x, float clip)
{
    // 输入: x 为真实浮点值，clip 为量化满量程。
    // 输出: Store 类型的饱和码值。
    // 用途: 执行裁剪、缩放、四舍五入和饱和，是 float -> qfloat 的核心步骤。
    if (clip <= 0.f) return static_cast<Store>(0);
    float xc = std::max(-clip, std::min(+clip, x));
    int   c  = static_cast<int>(std::lrint(xc * (static_cast<float>(Q()) / clip)));
    return static_cast<Store>(sat(c));
}

// Batch helpers (Matrix interop)
template<int NBITS>
matrix::Matrix< qfloat<NBITS> > quantize_matrix_to_qfloat(const matrix::Matrix<float>& in,
                                                          float clip)
{
    // 输入: in 为 float 矩阵，clip 为量化满量程。
    // 输出: 元素类型变为 qfloat<NBITS> 的新矩阵。
    // 用途: 批量准备量化输入数据。
    matrix::Matrix< qfloat<NBITS> > out(in.rows(), in.cols());
    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = qfloat<NBITS>::from_float(in[r][c], clip);
    return out;
}

template<int NBITS>
matrix::Matrix<float> dequantize_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in,
                                                     float clip_unused)
{
    // 输入: in 为 qfloat 矩阵，clip_unused 未实际使用。
    // 输出: 反量化后的 float 矩阵。
    // 用途: 批量查看或输出真实幅度近似值。
    (void)clip_unused;
    matrix::Matrix<float> out(in.rows(), in.cols());
    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = in[r][c].to_float();
    return out;
}

template<int NBITS>
matrix::Matrix<float> cast_matrix_from_qfloat(const matrix::Matrix< qfloat<NBITS> >& in)
{
    // 输入: in 为 qfloat 矩阵。
    // 输出: 每个元素等于对应 qfloat.code() 的 float 矩阵。
    // 用途: 调试或在“码值域”算法中观察内部整数表示。
    matrix::Matrix<float> out(in.rows(), in.cols());
    for (size_t r = 0; r < in.rows(); ++r)
        for (size_t c = 0; c < in.cols(); ++c)
            out[r][c] = static_cast<float>(static_cast<int>(in[r][c]));
    return out;
}

template<typename InputIt>
float compute_clip_from_ratio(InputIt first, InputIt last, float ratio)
{
    // 输入: [first, last) 为输入范围，ratio 为保留比例，范围会被夹到 [0,1]。
    // 输出: 按绝对值从大到小排序后，第 ceil(ratio * N) 个样本的幅度。
    // 用途: 用样本统计量自动选择 clip，减少过多饱和或分辨率浪费。
    std::vector<float> mags;
    for (auto it = first; it != last; ++it) {
        mags.push_back(std::fabs(static_cast<float>(*it)));
    }

    if (mags.empty()) return 0.0f;

    ratio = std::clamp(ratio, 0.0f, 1.0f);
    const std::size_t count = mags.size();
    std::size_t target = static_cast<std::size_t>(std::ceil(ratio * count));
    if (target == 0) target = 1;
    const std::size_t index = target - 1;

    std::nth_element(mags.begin(), mags.begin() + index, mags.end(), std::greater<float>());
    return mags[index];
}

} // namespace qfloat
