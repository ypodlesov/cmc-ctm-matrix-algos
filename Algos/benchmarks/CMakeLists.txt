cmake_minimum_required(VERSION 3.22)

project(benchmarks CXX)
set(CMAKE_CXX_FLAGS "-pedantic")

enable_testing()

include_directories("${PROJECT_SOURCE_DIR}/../src")

find_package(benchmark CONFIG REQUIRED)
include_directories(${benchmark_INCLUDE_DIRS})

add_executable(
    ${PROJECT_NAME}
    # QR_decomposition.cpp
    multiply.cpp
)

target_link_libraries(${PROJECT_NAME} benchmark::benchmark_main)
target_link_libraries(${PROJECT_NAME} ${QR_LIBRARY})

