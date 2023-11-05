cmake_minimum_required(VERSION 3.22)

project(fortran_examples)
enable_language(Fortran)


# if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#     set(dialect "-ffree-form -std=f2008 -fimplicit-none")
#     set(bounds "-fbounds-check")
# endif()
# if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#     set(dialect "-stand f08 -free -implicitnone")
#     set(bounds "-check bounds")
# endif()
# if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
#     set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
#     set(bounds "-C")
# endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

# file(GLOB_RECURSE modules src/matrix.f90)
# file(GLOB_RECURSE sources src/main.f90)
add_library(matrix STATIC src/matrix.f90)
add_library(lu STATIC src/lu.f90)
add_executable(main src/main.f90)
#the name of the external library located in /home/path/lib/libdir/lib is    libexternlib.so
target_link_libraries(main matrix) 
target_link_libraries(main lu) 
target_link_libraries(main blas) 