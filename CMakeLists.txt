cmake_minimum_required(VERSION 3.11)

project(QR_test CXX)
set(CMAKE_CXX_STANDARD 23)
enable_testing()

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

find_package(Boost 1.82.0 REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package(GTest REQUIRED)
if(GTest_FOUND)
  include_directories(BEFORE ${GTEST_INCLUDE_DIRS})
endif()

add_executable(
    QR_test 
    src/QR_test.cpp 
    src/QR.cpp
)

target_link_libraries(
    QR_test 
    gtest
    gtest_main
    rt
    pthread
)

