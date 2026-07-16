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
set(NEWCODE_INCLUDE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rx
)

# pipeline_runner 作为公共实现编入 newcode_frontend（apps 无需重复编译）

# ==============================
# 源文件自动收集（更耐用）
# ==============================
# 核心/通道/提取模块
file(GLOB CORE_SRCS    CONFIGURE_DEPENDS "src/common/core/*.cpp")
file(GLOB CHANNEL_SRCS CONFIGURE_DEPENDS "src/channel/*.cpp")
file(GLOB EXTRACT_SRCS CONFIGURE_DEPENDS "src/common/extract/*.cpp")
file(GLOB RX_BER_SRCS  CONFIGURE_DEPENDS "src/rx/ber/*.cpp")
file(GLOB COMMON_INTERLEAVER_SRCS CONFIGURE_DEPENDS "src/common/interleaver/*.cpp")
file(GLOB COMMON_PIPELINE_SRCS    CONFIGURE_DEPENDS "src/common/pipeline/*.cpp")
file(GLOB COMMON_UTILS_SRCS       CONFIGURE_DEPENDS "src/common/utils/*.cpp")
file(GLOB RX_DECODER_REGISTRY_SRCS CONFIGURE_DEPENDS "src/rx/decoder/registry/*.cpp")

# qfloat implementation sources (显式加入以生成模板显式实例化符号)
file(GLOB QFLOAT_SRCS CONFIGURE_DEPENDS "src/common/qfloat/*.cpp")

file(GLOB TX_BITGEN_SRCS CONFIGURE_DEPENDS "src/tx/bitgen/*.cpp")
file(GLOB TX_OFEC_SRCS   CONFIGURE_DEPENDS "src/tx/ofec/*.cpp")
file(GLOB TX_QAM_SRCS    CONFIGURE_DEPENDS "src/tx/qam/*.cpp")

file(GLOB RX_QAM_SRCS     CONFIGURE_DEPENDS "src/rx/qam/*.cpp")
file(GLOB RX_LLR_SRCS     CONFIGURE_DEPENDS "src/rx/llr/*.cpp")
file(GLOB RX_EXTRACT_SRCS CONFIGURE_DEPENDS "src/rx/extract/*.cpp")

# 两个解码变体
file(GLOB DECODER_PLAIN_SRCS   CONFIGURE_DEPENDS "src/rx/ofec/plain/*.cpp")
file(GLOB DECODER_EBCHPF_SRCS  CONFIGURE_DEPENDS "src/rx/ofec/ebchPF/*.cpp")

# 若存在公共译码源（两变体共享），自动并入两边
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/rx/ofec/common")
  file(GLOB DECODER_COMMON_SRCS CONFIGURE_DEPENDS "src/rx/ofec/common/*.cpp")
  list(APPEND DECODER_PLAIN_SRCS  ${DECODER_COMMON_SRCS})
  list(APPEND DECODER_EBCHPF_SRCS ${DECODER_COMMON_SRCS})
endif()

# 共享译码核心（plain / ebchPF 共用）
set(DECODER_SHARED_SRCS
  src/rx/ofec/ofec_row_decoder_core.cpp
  src/rx/ofec/ofec_tile_process.cpp
  src/rx/ofec/ofec_window_process.cpp
  src/rx/ofec/ofec_frame_decode.cpp
  src/rx/ofec/ofec_hard_decode_bch.cpp
)

# ==============================
# Shared "Frontend" Library（解码前/公共流水线）
# ==============================
add_library(newcode_frontend STATIC
  ${CORE_SRCS}
  ${CHANNEL_SRCS}
  ${EXTRACT_SRCS}
  ${COMMON_PIPELINE_SRCS}
  ${COMMON_INTERLEAVER_SRCS}
  ${COMMON_UTILS_SRCS}
  ${RX_DECODER_REGISTRY_SRCS}
  ${QFLOAT_SRCS}
  ${TX_BITGEN_SRCS}
  ${TX_OFEC_SRCS}
  ${TX_QAM_SRCS}
  ${RX_QAM_SRCS}
  ${RX_LLR_SRCS}
  ${RX_EXTRACT_SRCS}
  ${RX_BER_SRCS}
  src/rx/ofec/mux/mux_state_builder.cpp
  src/rx/ofec/mux/mux_siso_budget.cpp
  src/rx/ofec/mux/mux_config_validate.cpp
  src/rx/ofec/mux/mux_group_budget.cpp
  src/rx/ofec/mux/mux_group_config_validate.cpp
  src/rx/ofec/mux/mux_topology.cpp
  src/rx/ofec/mux/mux_matching.cpp
  src/rx/ofec/mux/mux_scheme_b_reconfig.cpp
  src/rx/ofec/mux/mux_scheme_c_staged.cpp
  src/rx/ofec/mux/mux_state_schedule_apply.cpp
)
target_include_directories(newcode_frontend PUBLIC ${NEWCODE_INCLUDE})
target_include_directories(newcode_frontend PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/ofec_sweep)
target_link_libraries(newcode_frontend PUBLIC Threads::Threads)

# ------- Compatibility fallbacks（临时兜底，便于平滑迁移）-------
# 编码侧依赖的 BCH：让 frontend 也能提供 parity/encode
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/common/bch")
  file(GLOB COMMON_BCH_SRCS CONFIGURE_DEPENDS "src/common/bch/*.cpp")
  target_sources(newcode_frontend PRIVATE ${COMMON_BCH_SRCS})
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/fec/encoder/bch_255_239.cpp")
  target_sources(newcode_frontend PRIVATE src/fec/encoder/bch_255_239.cpp)
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/bch_255_239.cpp")
  target_sources(newcode_frontend PRIVATE src/bch_255_239.cpp)
endif()

# 如果 ofec_encoder 仍在旧路径，也一起编进 frontend
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/ofec_encoder.cpp")
  target_sources(newcode_frontend PRIVATE src/ofec_encoder.cpp)
endif()

# ==============================
# Decoder Variants（解码后两套实现）
# ==============================
add_library(decoder_plain STATIC
  ${DECODER_PLAIN_SRCS}
  ${DECODER_SHARED_SRCS}
)
target_include_directories(decoder_plain PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(decoder_plain PUBLIC newcode_frontend)
target_sources(decoder_plain PRIVATE src/rx/decoder/factories/plain_decoder_factory.cpp)
target_sources(decoder_plain PRIVATE src/rx/decoder/factories/topk_pruned_decoder_factory.cpp)
target_sources(decoder_plain PRIVATE src/rx/decoder/factories/global_pair_decoder_factory.cpp)
target_sources(decoder_plain PRIVATE src/rx/decoder/factories/group_minima_decoder_factory.cpp)

add_library(decoder_ebchPF STATIC
  ${DECODER_EBCHPF_SRCS}
  ${DECODER_SHARED_SRCS}
)
target_include_directories(decoder_ebchPF PUBLIC ${NEWCODE_INCLUDE})
target_link_libraries(decoder_ebchPF PUBLIC newcode_frontend)
target_sources(decoder_ebchPF PRIVATE src/rx/decoder/factories/ebchPF_decoder_factory.cpp)

# 如果专门的模板实例化文件存在（建议放在 src/fec/decoder/common/）
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/rx/ofec/common/ofec_decoder_inst.cpp")
  target_sources(decoder_plain  PRIVATE src/rx/ofec/common/ofec_decoder_inst.cpp)
  target_sources(decoder_ebchPF PRIVATE src/rx/ofec/common/ofec_decoder_inst.cpp)
endif()

# ==============================
# Apps
# ==============================
# 保留现有两个入口，默认走 plain 变体
add_executable(ofec_single
  apps/ofec_single.cpp
  src/ofec_single/ofec_single_dual_writer.cpp
  src/ofec_single/ofec_single_params.cpp
  src/ofec_single/ofec_single_summary.cpp
  src/ofec_single/ofec_single_runner.cpp
)
target_link_libraries(ofec_single PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_sweep
  apps/ofec_sweep.cpp
  src/ofec_sweep/ofec_sweep_io.cpp
  src/ofec_sweep/ofec_sweep_scenarios.cpp
  src/ofec_sweep/ofec_sweep_runner.cpp
)
target_include_directories(ofec_sweep PRIVATE src/ofec_sweep)
target_link_libraries(ofec_sweep  PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_sweep2
  apps/ofec_sweep2.cpp
  src/ofec_sweep/ofec_sweep_io.cpp
  src/ofec_sweep/ofec_sweep_scenarios.cpp
  src/ofec_sweep/ofec_sweep_runner.cpp
)
target_include_directories(ofec_sweep2 PRIVATE src/ofec_sweep)
target_link_libraries(ofec_sweep2 PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_sweep3
  apps/ofec_sweep3.cpp
  src/ofec_sweep/ofec_sweep_io.cpp
  src/ofec_sweep/ofec_sweep_scenarios.cpp
  src/ofec_sweep/ofec_sweep_runner.cpp
)
target_include_directories(ofec_sweep3 PRIVATE src/ofec_sweep)
target_link_libraries(ofec_sweep3 PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_ber_window_probe
  apps/ofec_ber_window_probe.cpp
  src/ofec_single/ofec_single_dual_writer.cpp
  src/ofec_single/ofec_single_params.cpp
)
target_link_libraries(ofec_ber_window_probe PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_hybrid_regression_check
  apps/ofec_hybrid_regression_check.cpp
)
target_link_libraries(ofec_hybrid_regression_check PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_level56_shared_regression_check
  apps/ofec_level56_shared_regression_check.cpp
)
target_link_libraries(ofec_level56_shared_regression_check PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_fast_classifier_compare
  apps/ofec_fast_classifier_compare.cpp
)
target_link_libraries(ofec_fast_classifier_compare PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_two_stream_shared_decoder
  apps/ofec_two_stream_shared_decoder.cpp
  src/ofec_single/ofec_single_dual_writer.cpp
  src/ofec_single/ofec_single_params.cpp
  src/two_stream_shared/shared_bch_hard_decode64.cpp
  src/two_stream_shared/two_stream_shared_runner.cpp
)
target_link_libraries(ofec_two_stream_shared_decoder PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_two_stream_shared_sweep
  apps/ofec_two_stream_shared_sweep.cpp
  src/ofec_single/ofec_single_dual_writer.cpp
  src/ofec_single/ofec_single_params.cpp
  src/ofec_sweep/ofec_sweep_io.cpp
  src/ofec_sweep/ofec_sweep_scenarios.cpp
  src/ofec_sweep/ofec_sweep_runner.cpp
  src/two_stream_shared/shared_bch_hard_decode64.cpp
  src/two_stream_shared/two_stream_shared_runner.cpp
)
target_include_directories(ofec_two_stream_shared_sweep PRIVATE src/ofec_sweep)
target_link_libraries(ofec_two_stream_shared_sweep PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(ofec_two_stream_ber_window_probe
  apps/ofec_two_stream_ber_window_probe.cpp
  src/ofec_single/ofec_single_dual_writer.cpp
  src/ofec_single/ofec_single_params.cpp
  src/two_stream_shared/shared_bch_hard_decode64.cpp
  src/two_stream_shared/two_stream_shared_runner.cpp
)
target_link_libraries(ofec_two_stream_ber_window_probe PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(tpc_single
  apps/TPC_single.cpp
  src/tpc/tpc_decoder.cpp
  src/tpc/tpc_encoder.cpp
  src/tpc/tpc_pipline.cpp
  src/tpc/tpc_single_runner.cpp
)
target_link_libraries(tpc_single PRIVATE newcode_frontend decoder_plain decoder_ebchPF)

add_executable(tpc_sweep
  apps/tpc_sweep.cpp
  src/tpc/tpc_decoder.cpp
  src/tpc/tpc_encoder.cpp
  src/tpc/tpc_pipline.cpp
  src/tpc_sweep/tpc_sweep_io.cpp
  src/tpc_sweep/tpc_sweep_scenarios.cpp
  src/tpc_sweep/tpc_sweep_runner.cpp
)
target_link_libraries(tpc_sweep PRIVATE newcode_frontend decoder_plain decoder_ebchPF)
