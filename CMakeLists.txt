cmake_minimum_required(VERSION 3.13)
project(SomaticVariantCalling)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

# no warnings
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra -Wfatal-errors -msse4")

find_package(OpenMP REQUIRED)

include_directories(.)
include_directories(third_party)
include_directories(third_party/armadillo-10.3.0/include)
include_directories(third_party/spdlog/include)

set(CMAKE_CXX_FLAGS_DEBUG "-g -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#Progress Bar library
add_library(progress_bar third_party/progress_bar/progress_bar.cpp)

# Google Flags Library
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/gflags EXCLUDE_FROM_ALL)

# added later because glags has errors on GCC
string(APPEND CMAKE_CXX_FLAGS " -Werror")

# Armadillo Linear algebra library
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/armadillo-10.3.0 EXCLUDE_FROM_ALL)
target_compile_options(armadillo PUBLIC -Wno-unused-parameter)

# Spdlog logging library
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/spdlog EXCLUDE_FROM_ALL)

file(GLOB util_files util.cpp pileup_reader.cpp logger.cpp)
add_library(lib_util ${util_files})
target_link_libraries(lib_util spdlog::spdlog)

file(GLOB expectation_maximization_files expectation_maximization.cpp)
add_library(lib_expectation_maximization ${expectation_maximization_files})
target_link_libraries(lib_expectation_maximization lib_util)

file(GLOB similarity_matrix_files similarity_matrix.cpp)
add_library(lib_similarity_matrix ${similarity_matrix_files})
target_link_libraries(lib_similarity_matrix progress_bar lib_util OpenMP::OpenMP_CXX)

file(GLOB spectral_clustering_files spectral_clustering.cpp)
add_library(lib_spectral_clustering ${spectral_clustering_files})
target_link_libraries(lib_spectral_clustering armadillo lib_util  OpenMP::OpenMP_CXX)

file(GLOB variant_calling_files variant_calling.cpp variant_calling.cpp)
add_library(lib_variant_calling ${variant_calling_files})
target_link_libraries(lib_variant_calling lib_util)

add_executable(expectation_maximization expectation_maximization_main.cpp)
target_link_libraries(expectation_maximization lib_expectation_maximization lib_util gflags OpenMP::OpenMP_CXX)

add_executable(svc svc_main.cpp)
target_link_libraries(svc lib_expectation_maximization lib_variant_calling lib_similarity_matrix
    lib_spectral_clustering lib_util gflags spdlog::spdlog)

# TESTS
enable_testing()
include(GoogleTest)

add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/googletest EXCLUDE_FROM_ALL)
target_compile_options(gtest_main PRIVATE -w)
target_compile_options(gtest PRIVATE -w)

file(GLOB test_files "tests/*.cpp")

add_executable(tests ${test_files})
target_link_libraries(tests gtest_main gtest gmock lib_expectation_maximization lib_spectral_clustering)
target_include_directories(tests PRIVATE "include")

gtest_discover_tests(tests)

# copy the data into the binary directory, so that tests can reference it
file(COPY tests/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
