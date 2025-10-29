cmake_minimum_required(VERSION 3.16)
project(newcode LANGUAGES CXX)

# ==============================
# Toolchain & Options
# ==============================
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(Threads REQUIRED)

# 统一头文件路径
set(NEWCODE_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# ==============================
# Shared "Frontend" Library
# （解码前/公共流水线）
# ==============================
# 如果你尚未把这些文件迁到 src/frontend/，请把路径改回当前实际路径即可
add_library(newcode_frontend STATIC
  src/frontend/bitgen.cpp
  src/frontend/ofec_encoder.cpp
  src/frontend/qam.cpp
  src/frontend/awgn.cpp
  src/frontend/qam_llr.cpp
  src/frontend/ofec_llr_matrix.cpp
  src/frontend/info_extract.cpp
  src/frontend/llr_qpack.cpp
  src/frontend/llr_known_prefix.cpp
  src/frontend/qfloat.cpp
  src/frontend/interleaver_identity.cpp
  src/frontend/interleaver_ofec.cpp
  # 如果你希望把 pipeline_runner 作为公共组件复用，请取消下一行注释并把文件移动到 frontend
  # src/frontend/pipeline_runner.cpp
)
target_include_directories(newcode_frontend PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(newcode_frontend PUBLIC Threads::Threads)

# ==============================
# Decoder Variants
# （解码后两套实现）
# ==============================
# —— Plain 变体：无交织/或你的“非 ebchPF”路径 —— #
add_library(decoder_plain STATIC
  src/decoders/plain/ofec_decoder.cpp
  src/decoders/plain/decoder_core.cpp
  src/decoders/plain/ofec_decoder_hard.cpp
  src/decoders/plain/chase256.cpp
  src/decoders/plain/llr_to_bit.cpp
  src/decoders/plain/bch_255_239.cpp
  src/decoders/plain/decoder_factory.cpp   # 你的工厂（make_decoder("plain")）
  # 如有其它专属实现文件，在此追加
)
target_include_directories(decoder_plain PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(decoder_plain PUBLIC newcode_frontend)

# —— ebchPF 变体：有交织/你的 ebchPF 路径 —— #
add_library(decoder_ebchPF STATIC
  src/decoders/ebchPF/ofec_decoder.cpp
  src/decoders/ebchPF/decoder_core.cpp
  src/decoders/ebchPF/ofec_decoder_hard.cpp
  src/decoders/ebchPF/chase256.cpp
  src/decoders/ebchPF/llr_to_bit.cpp
  src/decoders/ebchPF/bch_255_239.cpp
  src/decoders/ebchPF/decoder_factory.cpp  # 你的工厂（make_decoder("ebchPF")）
  # 如有其它专属实现文件，在此追加
)
target_include_directories(decoder_ebchPF PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(decoder_ebchPF PUBLIC newcode_frontend)

# ==============================
# Apps
# ==============================
# 方案A：两个新的 demo 入口（建议）
add_executable(ofec_app_plain
  apps/ofec_app_plain.cpp
  # 若 pipeline_runner 仍在 src/ 下而非 frontend，请把它编进可执行：
  src/pipeline_runner.cpp
)
target_link_libraries(ofec_app_plain
  PRIVATE newcode_frontend decoder_plain
)

add_executable(ofec_app_interleaved
  apps/ofec_app_interleaved.cpp
  # 同上，如未迁移：
  src/pipeline_runner.cpp
)
target_link_libraries(ofec_app_interleaved
  PRIVATE newcode_frontend decoder_ebchPF
)

# 方案B：保留你现有的两个入口，默认走 plain 变体（你可按需改为 ebchPF）
add_executable(ofec_single apps/ofec_single.cpp src/pipeline_runner.cpp)
target_link_libraries(ofec_single PRIVATE newcode_frontend decoder_plain)

add_executable(ofec_sweep  apps/ofec_sweep.cpp  src/pipeline_runner.cpp)
target_link_libraries(ofec_sweep  PRIVATE newcode_frontend decoder_plain)

# ==============================
# Tests
# ==============================
include(CTest)
if(BUILD_TESTING)
  # 这些测试大多只依赖“前处理”与基础类型即可
  add_executable(test_matrix tests/test_matrix.cpp)
  target_link_libraries(test_matrix PRIVATE newcode_frontend)
  add_test(NAME MatrixTest COMMAND test_matrix)

  add_executable(test_qam tests/test_qam.cpp)
  target_link_libraries(test_qam PRIVATE newcode_frontend)
  add_test(NAME QAMTest COMMAND test_qam)

  add_executable(test_awgn tests/test_awgn.cpp)
  target_link_libraries(test_awgn PRIVATE newcode_frontend)
  add_test(NAME AWGNTest COMMAND test_awgn)

  add_executable(test_qam_llr tests/test_qam_llr.cpp)
  target_link_libraries(test_qam_llr PRIVATE newcode_frontend)
  add_test(NAME QAMLLRTest COMMAND test_qam_llr)

  add_executable(llr_reshape_demo tests/llr_reshape_demo.cpp)
  target_link_libraries(llr_reshape_demo PRIVATE newcode_frontend)
  add_test(NAME LLRReshapeTest COMMAND llr_reshape_demo)

  add_executable(test_interleaver tests/test_interleaver.cpp)
  target_link_libraries(test_interleaver PRIVATE newcode_frontend)
  add_test(NAME InterleaverTest COMMAND test_interleaver)

  # BCH 专项测试通常依赖译码实现，这里链接到一个变体（任选其一）
  add_executable(test_bch tests/test_bch.cpp)
  target_link_libraries(test_bch PRIVATE newcode_frontend decoder_plain)
  add_test(NAME BCHTest COMMAND test_bch)
endif()
