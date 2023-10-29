cmake_minimum_required(VERSION 3.22)

project(algos CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

set(QR_LIBRARY QR_decomposition_lib)
add_subdirectory(src)

if (DEFINED test)
add_subdirectory(tests)
endif()

if (DEFINED benchmark)
add_subdirectory(benchmarks)
endif()
