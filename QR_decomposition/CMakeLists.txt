cmake_minimum_required(VERSION 3.22)

project(QR_decomposition CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

set(QR_LIBRARY QR_decomposition_lib)
add_subdirectory(src)

if (DEFINED test)
add_subdirectory(test)
endif()

if (DEFINED benchmark)
add_subdirectory(benchmark)
endif()
