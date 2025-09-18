cmake_minimum_required(VERSION 3.16)
project(newcode LANGUAGES CXX)

# 使用 C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 警告选项
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 可执行程序，包含 main.cpp 和 newcode.cpp
add_executable(newcode_app
  src/main.cpp
  src/bitgen.cpp
)

# 头文件目录
target_include_directories(newcode_app PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
