cmake_minimum_required(VERSION 3.22)

project(benchmarks CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

enable_testing()

include_directories("${PROJECT_SOURCE_DIR}/../src")

find_package(benchmark CONFIG REQUIRED)
include_directories(${benchmark_INCLUDE_DIRS})

find_package(GTest REQUIRED)
include_directories(${GTest_INCLUDE_DIRS})

add_executable(
    ${PROJECT_NAME}
    qr_decomposition.cpp
    # multiply.cpp
)

target_link_libraries(${PROJECT_NAME} benchmark::benchmark)
target_link_libraries(${PROJECT_NAME} benchmark::benchmark_main)
target_link_libraries(${PROJECT_NAME} ${QR_LIBRARY})

