cmake_minimum_required(VERSION 3.22)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

# enable_testing()

include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/src")

# find_package(GTest REQUIRED)
# include_directories(${GTest_INCLUDE_DIRS})

add_executable(
    ${PROJECT_NAME}
    main.cpp
)
    
# target_link_libraries(${PROJECT_NAME} GTest::GTest)
# target_link_libraries(${PROJECT_NAME} GTest::gtest_main)
# target_link_libraries(${PROJECT_NAME} GTest::gmock)
# target_link_libraries(${PROJECT_NAME} GTest::gmock_main)
target_link_libraries(${PROJECT_NAME} ${QR_LIBRARY})