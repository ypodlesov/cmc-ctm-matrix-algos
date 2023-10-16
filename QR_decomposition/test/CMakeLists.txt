cmake_minimum_required(VERSION 3.22)

project(QR_decomposition_test CXX)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

enable_testing()

include_directories("${PROJECT_SOURCE_DIR}/../src")

find_package(GTest CONFIG REQUIRED)
include_directories(${GTest_INCLUDE_DIRS})

add_executable(
    ${PROJECT_NAME} 
    QR_decomposition.cpp
)
    
target_link_libraries(${PROJECT_NAME} gtest::gtest)
target_link_libraries(${PROJECT_NAME} ${QR_LIBRARY})