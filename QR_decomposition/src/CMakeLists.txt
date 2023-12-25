cmake_minimum_required(VERSION 3.22)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)


include_directories("${PROJECT_SOURCE_DIR}/include")
add_library(${QR_LIBRARY} SHARED)

set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
target_link_libraries(${QR_LIBRARY} PUBLIC -lopenblas)
target_compile_options(${QR_LIBRARY} PRIVATE -pedantic -march=native)

target_sources(${QR_LIBRARY}
    PUBLIC
        matrix.cpp
        vector.cpp
        qr_decomposition.cpp
)
