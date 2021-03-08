cmake_minimum_required(VERSION 3.13)
project(SomaticVariantCalling)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

find_package(OpenMP REQUIRED)

include_directories(.)
include_directories(third_party)

set(CMAKE_CXX_FLAGS_DEBUG "-g -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

#Progress Bar library
add_library(progress_bar third_party/progress_bar/progress_bar.cpp)

# Google Flags Library
add_subdirectory(${PROJECT_SOURCE_DIR}/third_party/gflags EXCLUDE_FROM_ALL)

file(GLOB similarity_matrix_files util.cpp)
add_library(lib_similarity_matrix ${similarity_matrix_files})

add_executable(similarity_matrix similarity_matrix_main.cpp)
target_link_libraries(similarity_matrix lib_similarity_matrix gflags OpenMP::OpenMP_CXX progress_bar)