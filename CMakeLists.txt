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
# 核心库：收纳所有通用实现源码，供主程�?测试复用
# ------------------------------------------------------------------------------
add_library(newcode_core STATIC
  src/bitgen.cpp
  # NOTE: OVERALLLEN is currently unused in ofec_encoder.cpp; mark [[maybe_unused]] if intentional.
  src/ofec_encoder.cpp
  src/bch_255_239.cpp
  src/qam.cpp
  src/awgn.cpp
  src/qam_llr.cpp
  src/ofec_llr_matrix.cpp
  src/ofec_decoder.cpp
  src/info_extract.cpp
  src/ber.cpp
  src/llr_to_bit.cpp
  src/llr_qpack.cpp
  src/llr_known_prefix.cpp
  src/chase256.cpp         # �?必须加入：提�?chase_decode_256<T>() 的定义与显式实例�?  # 如果你把 qfloat 做成 .cpp（非 header-only），也需要在这里加：
    src/pipeline_runner.cpp
)
target_include_directories(newcode_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# ------------------------------------------------------------------------------
# 主程�?# ------------------------------------------------------------------------------
add_executable(newcode_app
  src/main.cpp
)
target_link_libraries(newcode_app PRIVATE newcode_core)

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
  add_executable(test_matrix tests/test_matrix.cpp)
  target_link_libraries(test_matrix PRIVATE newcode_core)
  add_test(NAME MatrixTest COMMAND test_matrix)

  add_executable(test_bch tests/test_bch.cpp)
  target_link_libraries(test_bch PRIVATE newcode_core)
  add_test(NAME BCHTest COMMAND test_bch)

  add_executable(test_qam tests/test_qam.cpp)
  target_link_libraries(test_qam PRIVATE newcode_core)
  add_test(NAME QAMTest COMMAND test_qam)

  add_executable(test_awgn tests/test_awgn.cpp)
  target_link_libraries(test_awgn PRIVATE newcode_core)
  add_test(NAME AWGNTest COMMAND test_awgn)

  add_executable(test_qam_llr tests/test_qam_llr.cpp)
  target_link_libraries(test_qam_llr PRIVATE newcode_core)
  add_test(NAME QAMLLRTest COMMAND test_qam_llr)

  add_executable(llr_reshape_demo tests/llr_reshape_demo.cpp)
  target_link_libraries(llr_reshape_demo PRIVATE newcode_core)
  add_test(NAME LLRReshapeTest COMMAND llr_reshape_demo)
endif()
