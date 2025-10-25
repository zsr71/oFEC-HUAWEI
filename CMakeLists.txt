cmake_minimum_required(VERSION 3.16)
project(newcode LANGUAGES CXX)

# C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Warnings
if(MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# pthread / std::thread
find_package(Threads REQUIRED)

# ------------------------------------------------------------------------------
# Core library：收纳所有算法/流水线实现
# ------------------------------------------------------------------------------
add_library(newcode_core STATIC
  src/bitgen.cpp
  src/ofec_encoder.cpp
  src/bch_255_239.cpp
  src/qam.cpp
  src/awgn.cpp
  src/qam_llr.cpp
  src/ofec_llr_matrix.cpp
  src/ofec_decoder.cpp
  src/decoder_core.cpp
  src/ofec_decoder_hard.cpp
  src/info_extract.cpp
  src/ber.cpp
  src/llr_to_bit.cpp
  src/llr_qpack.cpp
  src/llr_known_prefix.cpp
  src/chase256.cpp
  src/pipeline_runner.cpp
)
target_include_directories(newcode_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(newcode_core PUBLIC Threads::Threads)

# ------------------------------------------------------------------------------
# Apps：两个入口（单次调试 / 并行扫参）
# ------------------------------------------------------------------------------
add_executable(ofec_single apps/ofec_single.cpp)
target_link_libraries(ofec_single PRIVATE newcode_core)

add_executable(ofec_sweep  apps/ofec_sweep.cpp)
target_link_libraries(ofec_sweep  PRIVATE newcode_core)

# （已删除 newcode_app/src/main.cpp，避免多个 main）
# 如果以后还想保留旧入口，请把旧 main 改名并放到 apps/ 下，作为单独的可执行去构建。

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
