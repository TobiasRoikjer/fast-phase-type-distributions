cmake_minimum_required(VERSION 3.7)
project(amazephase C)

set(CMAKE_C_FLAGS "-std=c11 ${CMAKE_C_FLAGS}")
set(CMAKE_C_STANDARD 11)

set(CMAKE_BINARY_DIR .)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

IF(CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O3")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

SET(SOURCE_FILES src/bbst.h src/utils.h src/coal.h src/sampling.h src/dist.h src/bbst.c src/coal.c src/dist.c src/utils.c src/sampling.c)
add_library(amazephase ${SOURCE_FILES})
add_executable(run_lmvu ./cmd/run_lmvu.c)
add_executable(run_build_kingman ./cmd/run_build_kingman.c)
add_executable(run_cmp_states_rw ./cmd/run_cmp_states_rw.c)
add_executable(run_time ./cmd/run_time.c)
add_executable(run_cdf_consts ./cmd/run_cdf_consts.c)
add_executable(run_pdf_consts_mut ./cmd/run_pdf_consts_mut.c)
add_executable(run_pure_cutoff ./cmd/run_pure_cutoff.c)
add_executable(run_cutoff ./cmd/run_cutoff.c)
add_executable(run_left_prob ./cmd/run_left_prob.c)
add_executable(run_im_cdf ./cmd/run_im_cdf.c)
add_executable(run_im_cdf_full ./cmd/run_im_cdf_full.c)
add_executable(run_im_mat ./cmd/run_im_mat.c)
add_executable(run_im_exp ./cmd/run_im_exp.c)
add_executable(run_im_full_prob_exp ./cmd/run_im_full_prob_exp.c)
add_executable(run_im_full_prob_cdf ./cmd/run_im_full_prob_cdf.c)
add_executable(run_num_coals_prob ./cmd/run_num_coals_prob.c)
add_executable(run_im_ss_state_size ./cmd/run_im_ss_state_size.c)
add_executable(run_disc ./cmd/run_disc.c)
add_executable(test_bbst src/test/test_bbst.c)
add_executable(test_utils src/test/test_utils.c)
add_executable(test_coal src/test/test_coal.c)
add_executable(test_sampling src/test/test_sampling.c)

find_package(GSL REQUIRED)
target_link_libraries(amazephase m gsl gslcblas)

target_link_libraries(run_lmvu amazephase)
target_link_libraries(run_build_kingman amazephase)
target_link_libraries(run_cmp_states_rw amazephase)
target_link_libraries(run_time amazephase)
target_link_libraries(run_cdf_consts amazephase)
target_link_libraries(run_pdf_consts_mut amazephase)
target_link_libraries(run_pure_cutoff amazephase)
target_link_libraries(run_cutoff amazephase)
target_link_libraries(run_left_prob amazephase)
target_link_libraries(run_im_cdf amazephase)
target_link_libraries(run_im_cdf_full amazephase)
target_link_libraries(run_im_mat amazephase)
target_link_libraries(run_im_exp amazephase)
target_link_libraries(run_im_full_prob_exp amazephase)
target_link_libraries(run_im_full_prob_cdf amazephase)
target_link_libraries(run_num_coals_prob amazephase)
target_link_libraries(run_im_ss_state_size amazephase)
target_link_libraries(run_disc amazephase)
target_link_libraries(test_bbst amazephase)
target_link_libraries(test_utils amazephase)
target_link_libraries(test_coal amazephase)
target_link_libraries(test_sampling amazephase)
