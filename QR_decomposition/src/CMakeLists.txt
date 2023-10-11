cmake_minimum_required(VERSION 3.22)

project(QR_decomposition CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

add_library(
    SHARED
    ${PROJECT_NAME} 
    matrix.cpp
    vector.cpp
    QR_decomposition.cpp
)