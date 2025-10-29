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

# 小助手：按存在情况选 pipeline_runner 源文件位置
set(PIPELINE_RUNNER_SRC "")
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/pipeline_runner.cpp")
  set(PIPELINE_RUNNER_SRC src/pipeline_runner.cpp)
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/frontend/pipeline_runner.cpp")
  set(PIPELINE_RUNNER_SRC src/frontend/pipeline_runner.cpp)
endif()

# ==============================
# 源文件自动收集（更耐用）
# ==============================
# 解码前 / 公共流水线
file(GLOB FRONTEND_SRCS CONFIGURE_DEPENDS
  "src/frontend/*.cpp"
)

# 两个解码变体
file(GLOB DECODER_PLAIN_SRCS   CONFIGURE_DEPENDS "src/decoders/plain/*.cpp")
file(GLOB DECODER_EBCHPF_SRCS  CONFIGURE_DEPENDS "src/decoders/ebchPF/*.cpp")

# 若存在公共译码源（两变体共享），自动并入两边
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/decoder_common")
  file(GLOB DECODER_COMMON_SRCS CONFIGURE_DEPENDS "src/decoder_common/*.cpp")
  list(APPEND DECODER_PLAIN_SRCS  ${DECODER_COMMON_SRCS})
  list(APPEND DECODER_EBCHPF_SRCS ${DECODER_COMMON_SRCS})
endif()

# ==============================
# Shared "Frontend" Library（解码前/公共流水线）
# ==============================
add_library(newcode_frontend STATIC
  ${FRONTEND_SRCS}
)
target_include_directories(newcode_frontend PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(newcode_frontend PUBLIC Threads::Threads)

# ==============================
# Decoder Variants（解码后两套实现）
# ==============================
# —— Plain 变体：无交织/非 ebchPF —— #
add_library(decoder_plain STATIC
  ${DECODER_PLAIN_SRCS}
)
target_include_directories(decoder_plain PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(decoder_plain PUBLIC newcode_frontend)

# —— ebchPF 变体：有交织/你的 ebchPF 路径 —— #
add_library(decoder_ebchPF STATIC
  ${DECODER_EBCHPF_SRCS}
)
target_include_directories(decoder_ebchPF PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(decoder_ebchPF PUBLIC newcode_frontend)
# ------- Compatibility fallbacks: use legacy file locations if still present -------
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/bch_255_239.cpp")
  target_sources(decoder_plain  PRIVATE src/bch_255_239.cpp)
  target_sources(decoder_ebchPF PRIVATE src/bch_255_239.cpp)
endif()

# 如果 ofec_encoder 还没移动到 frontend，也先把旧路径编入（以便 llr_reshape_demo 链接过）
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/ofec_encoder.cpp")
  target_sources(newcode_frontend PRIVATE src/ofec_encoder.cpp)
endif()

# ==============================
# Apps
# ==============================
# 新的 demo 入口
add_executable(ofec_app_plain
  apps/ofec_app_plain.cpp
  ${PIPELINE_RUNNER_SRC}
)
target_link_libraries(ofec_app_plain
  PRIVATE newcode_frontend decoder_plain
)

add_executable(ofec_app_interleaved
  apps/ofec_app_interleaved.cpp
  ${PIPELINE_RUNNER_SRC}
)
target_link_libraries(ofec_app_interleaved
  PRIVATE newcode_frontend decoder_ebchPF
)

# 保留你现有的两个入口，默认走 plain 变体（可按需改）
add_executable(ofec_single apps/ofec_single.cpp ${PIPELINE_RUNNER_SRC})
target_link_libraries(ofec_single PRIVATE newcode_frontend decoder_plain)

add_executable(ofec_sweep  apps/ofec_sweep.cpp  ${PIPELINE_RUNNER_SRC})
target_link_libraries(ofec_sweep  PRIVATE newcode_frontend decoder_plain)

# ==============================
# Tests
# ==============================
include(CTest)
if(BUILD_TESTING)
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
