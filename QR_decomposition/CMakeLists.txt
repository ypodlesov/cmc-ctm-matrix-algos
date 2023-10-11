cmake_minimum_required(VERSION 3.22)

project(QR_decomposition CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

add_subdirectory(src)
add_subdirectory(test)

# find_package(benchmark CONFIG REQUIRED)
# include_directories(${benchmark_INCLUDE_DIRS})

# target_link_libraries(${PROJECT_NAME} gtest)
# target_link_libraries(${PROJECT_NAME} benchmark::benchmark_main)

