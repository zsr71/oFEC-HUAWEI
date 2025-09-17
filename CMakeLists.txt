cmake_minimum_required(VERSION 3.16)
project(newcode LANGUAGES CXX)

# C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 警告更严格（可选）
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 库
add_library(newcode STATIC
  src/newcode.cpp
)
target_include_directories(newcode PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# 可执行程序
add_executable(newcode_app src/main.cpp)
target_link_libraries(newcode_app PRIVATE newcode)
