#pragma once
#include <string>

namespace newcode {

// 一个简单的版本结构
struct Version { int major, minor, patch; };
Version version();

// 随便放个演示函数：计算整数奇偶校验位
int parity(unsigned x);

// 打个招呼
std::string hello(const std::string& name);

} // namespace newcode
