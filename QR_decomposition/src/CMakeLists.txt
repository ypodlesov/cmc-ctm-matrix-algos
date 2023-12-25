cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-pedantic")

add_library(${QR_LIBRARY} SHARED)

find_package(LAPACK REQUIRED)
target_link_libraries(${QR_LIBRARY} ${LAPACK_LIBRARIES})

target_sources(${QR_LIBRARY}
    PUBLIC
        matrix.cpp
        vector.cpp
        qr_decomposition.cpp
)