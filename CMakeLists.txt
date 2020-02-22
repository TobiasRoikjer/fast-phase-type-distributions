# cmake_minimum_required(VERSION <specify CMake version here>)
project(amazephase C)

set(CMAKE_C_FLAGS "-std=c11 ${CMAKE_C_FLAGS}")
set(CMAKE_C_STANDARD 11)

set(CMAKE_BINARY_DIR .)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(amazephase bbst.h utils.h bbst.c mat.c mat.h coal.c coal.h phdist.c phdist.h dist.c dist.h utils.c sampling.c sampling.h)
add_executable(run_lmvu run_lmvu.c bbst.h utils.h bbst.c mat.c mat.h coal.c coal.h phdist.c phdist.h dist.c dist.h utils.c sampling.c sampling.h)
add_executable(test_bbst test/test_bbst.c bbst.h utils.h bbst.c mat.c mat.h coal.c coal.h dist.c dist.h utils.c sampling.c sampling.h)
add_executable(test_utils test/test_utils.c utils.h utils.c sampling.c sampling.h)
add_executable(test_coal test/test_coal.c bbst.h utils.h mat.h coal.h phdist.h bbst.c mat.c coal.c phdist.c dist.c dist.h utils.c sampling.c sampling.h)
add_executable(test_dist test/test_sampling.c bbst.h utils.h mat.h coal.h phdist.h bbst.c mat.c coal.c phdist.c dist.c dist.h utils.c sampling.c sampling.h)
target_link_libraries(amazephase m)
target_link_libraries(test_dist m)
