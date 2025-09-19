cmake_minimum_required(VERSION 3.16)
project(newcode LANGUAGES CXX)

# C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Warnings
if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# ------------------------------------------------------------------------------
# 核心库：收纳所有通用实现源码，供主程序/测试复用
# ------------------------------------------------------------------------------
add_library(newcode_core STATIC
  src/bitgen.cpp
  src/ofec_encoder.cpp
  src/bch_255_239.cpp
  src/qam.cpp
  src/awgn.cpp
  src/qam_llr.cpp
)
target_include_directories(newcode_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ------------------------------------------------------------------------------
# 主程序
# ------------------------------------------------------------------------------
add_executable(newcode_app
  src/main.cpp
)
target_link_libraries(newcode_app PRIVATE newcode_core)

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------
include(CTest)         # 等价于 enable_testing() 且支持 BUILD_TESTING 开关
if(BUILD_TESTING)
  # Matrix test
  add_executable(test_matrix tests/test_matrix.cpp)
  target_link_libraries(test_matrix PRIVATE newcode_core)
  add_test(NAME MatrixTest COMMAND test_matrix)

  # BCH test
  add_executable(test_bch tests/test_bch.cpp)
  target_link_libraries(test_bch PRIVATE newcode_core)
  add_test(NAME BCHTest COMMAND test_bch)

  # QAM test
  add_executable(test_qam tests/test_qam.cpp)
  target_link_libraries(test_qam PRIVATE newcode_core)
  add_test(NAME QAMTest COMMAND test_qam)

  # AWGN test
  add_executable(test_awgn tests/test_awgn.cpp)
  target_link_libraries(test_awgn PRIVATE newcode_core)
  add_test(NAME AWGNTest COMMAND test_awgn)

  # 新测试
  add_executable(test_qam_llr tests/test_qam_llr.cpp)
  target_link_libraries(test_qam_llr PRIVATE newcode_core)
  add_test(NAME QAMLLRTest COMMAND test_qam_llr)

endif()
