cmake_minimum_required(VERSION 3.13)
project(SomaticVariantCalling)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

find_package(OpenMP REQUIRED)
# no warnings
string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra -Wfatal-errors -msse4 ${OpenMP_CXX_FLAGS}")

include_directories(.)
include_directories(third_party)
include_directories(third_party/armadillo-10.3.0/include)
include_directories(third_party/spdlog/include)

set(CMAKE_CXX_FLAGS_DEBUG "-g -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#Progress Bar library
add_library(progress_bar third_party/progress_bar/progress_bar.cpp)
target_compile_options(progress_bar PUBLIC -fPIC)

# Google Flags Library
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/gflags EXCLUDE_FROM_ALL)

# added later because glags has errors on GCC
string(APPEND CMAKE_CXX_FLAGS " -Werror")

# Armadillo Linear algebra library
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/armadillo-10.3.0 EXCLUDE_FROM_ALL)
target_compile_options(armadillo PUBLIC -Wno-unused-parameter)

# Spdlog logging library
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/spdlog EXCLUDE_FROM_ALL)

# BAMTools library for parsing BAM files
include_directories(${PROJECT_SOURCE_DIR}/third_party/bamtools/src)
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/bamtools EXCLUDE_FROM_ALL)


file(GLOB util_files util/*.cpp)
add_library(util ${util_files})
target_link_libraries(util spdlog::spdlog armadillo BamTools progress_bar OpenMP::OpenMP_CXX)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(util stdc++fs)
endif()

file(GLOB expectation_maximization_files expectation_maximization.cpp)
add_library(expectation_maximization ${expectation_maximization_files})
target_link_libraries(expectation_maximization util OpenMP::OpenMP_CXX)

file(GLOB similarity_matrix_files similarity_matrix.cpp)
add_library(similarity_matrix ${similarity_matrix_files})
target_link_libraries(similarity_matrix progress_bar util OpenMP::OpenMP_CXX)

file(GLOB spectral_clustering_files spectral_clustering.cpp)
add_library(spectral_clustering ${spectral_clustering_files})
target_link_libraries(spectral_clustering armadillo expectation_maximization similarity_matrix util  OpenMP::OpenMP_CXX)

file(GLOB variant_calling_files variant_calling.cpp variant_calling.cpp)
add_library(variant_calling ${variant_calling_files})
target_link_libraries(variant_calling util)

file(GLOB pileup_files pileup.cpp pileup.cpp)
add_library(pileup_lib ${pileup_files})
target_link_libraries(pileup_lib util)

add_executable(silver silver_main.cpp)
target_link_libraries(silver expectation_maximization variant_calling similarity_matrix
    spectral_clustering util gflags spdlog::spdlog BamTools)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(silver stdc++fs)
endif()

add_executable(pileup pileup_main.cpp)
target_link_libraries(pileup pileup_lib util gflags spdlog::spdlog BamTools)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(pileup stdc++fs)
endif()

# TESTS
enable_testing()
include(GoogleTest)

add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/googletest EXCLUDE_FROM_ALL)
target_compile_options(gtest_main PRIVATE -w)
target_compile_options(gtest PRIVATE -w)

file(GLOB test_files "tests/*.cpp")

add_executable(tests ${test_files})
target_link_libraries(tests gtest_main gtest gmock expectation_maximization pileup_lib spectral_clustering
    variant_calling)
target_include_directories(tests PRIVATE "include")

gtest_discover_tests(tests)

# copy the data into the binary directory, so that tests can reference it
file(COPY tests/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
