cmake_minimum_required(VERSION 3.16)
project(new_float_only LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(Threads REQUIRED)

set(NEW_FLOAT_ONLY_INCLUDE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common
  ${CMAKE_CURRENT_SOURCE_DIR}/src/rx
)

add_library(new_float_only_core STATIC
  src/common/bch/bch_255_239_decode.cpp
  src/common/bch/bch_255_239_encode.cpp
  src/common/bch/bch_255_239_parity.cpp
  src/common/utils/ensure_dir.cpp
  src/common/utils/hard_bits_to_llr_matrix.cpp
  src/common/utils/now_stamp.cpp
  src/common/utils/prepare_log_file.cpp
  src/channel/awgn.cpp
  src/channel/ebn0_to_sigma.cpp
  src/pipeline/pipeline_runner.cpp
  src/rx/ber/ber.cpp
  src/rx/extract/info_extract.cpp
  src/rx/llr/llr_known_prefix.cpp
  src/rx/llr/ofec_llr_matrix.cpp
  src/rx/ofec/ofec_frame_decode.cpp
  src/rx/ofec/ofec_hard_decode_bch.cpp
  src/rx/ofec/ofec_row_decoder_core.cpp
  src/rx/ofec/ofec_tile_process.cpp
  src/rx/ofec/ofec_window_process.cpp
  src/rx/ofec/plain/chase256_plain.cpp
  src/rx/qam/qam_llr.cpp
  src/tx/bitgen/bitgen.cpp
  src/tx/ofec/ofec_encoder.cpp
  src/tx/qam/qam.cpp
  src/ofec_single/ofec_single_runner.cpp
  src/ofec_single/ofec_single_dual_writer.cpp
  src/seed_sweep/ebn0_sweep_runner.cpp
  src/seed_sweep/seed_sweep_runner.cpp
  src/seed_sweep/sweep_task_runner.cpp
)

target_include_directories(new_float_only_core PUBLIC ${NEW_FLOAT_ONLY_INCLUDE})
target_link_libraries(new_float_only_core PUBLIC Threads::Threads)

add_executable(ofec_single_float apps/ofec_single_float.cpp)
target_link_libraries(ofec_single_float PRIVATE new_float_only_core)

add_executable(ofec_seed_sweep_float apps/ofec_seed_sweep_float.cpp)
target_link_libraries(ofec_seed_sweep_float PRIVATE new_float_only_core)

add_executable(ofec_alpha_beta_sweep_float apps/ofec_alpha_beta_sweep_float.cpp)
target_link_libraries(ofec_alpha_beta_sweep_float PRIVATE new_float_only_core)

add_executable(ofec_ebn0_sweep_float apps/ofec_ebn0_sweep_float.cpp)
target_link_libraries(ofec_ebn0_sweep_float PRIVATE new_float_only_core)
