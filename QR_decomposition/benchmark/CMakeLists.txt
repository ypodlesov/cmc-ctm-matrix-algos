cmake_minimum_required(VERSION 3.22)

project(QR_decomposition_benchmark CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")
# set(BENCHMARK_OUT_FORMAT csv)

enable_testing()

include_directories("${PROJECT_SOURCE_DIR}/..")
include_directories("${PROJECT_SOURCE_DIR}/../src")


find_package(benchmark CONFIG REQUIRED)
include_directories(${benchmark_INCLUDE_DIRS})

add_executable(
    ${PROJECT_NAME} 
    QR_decomposition.cpp
    ../src/matrix.cpp
    ../src/vector.cpp
    ../src/QR_decomposition.cpp
)

target_link_libraries(${PROJECT_NAME} benchmark::benchmark_main)

