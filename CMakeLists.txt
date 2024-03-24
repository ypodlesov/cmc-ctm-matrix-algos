cmake_minimum_required(VERSION 3.14)
project(shad-cpp)

set(CMAKE_CXX_STANDARD             20)
set(CMAKE_MODULE_PATH              "${CMAKE_SOURCE_DIR}/tools/cmake")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
set(CMAKE_EXPORT_COMPILE_COMMANDS  ON)

find_package(Catch REQUIRED)
find_package(TBB REQUIRED)

find_package(PNG)
find_package(JPEG)
find_package(FFTW)

set(Clang_DIR "/usr/lib/llvm-16/lib/cmake/clang")
set(LLVM_DIR "/usr/lib/llvm-16/lib/cmake/llvm")
find_package(Clang CONFIG)

include_directories(tools/commons)
include(tools/cmake/TestSolution.cmake)
include(tools/cmake/BuildFlags.cmake)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/private)
  add_subdirectory(private)
endif()

function(add_if_exists name)
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name})
    add_subdirectory(${name})
  endif()
endfunction()

add_if_exists(qr_decomposition)

# if (PNG_FOUND AND JPEG_FOUND AND FFTW_FOUND)
#   add_if_exists(jpeg-decoder)
# elseif()
#   message(STATUS "jpeg-decoder disabled. PNG:${PNG_FOUND} JPEG:${JPEG_FOUND} FFTW:${FFTW_FOUND}")
# endif()

# add_if_exists(spinlock)
# add_if_exists(rw-spinlock)
# add_if_exists(rw-counter)
# add_if_exists(fast-queue)
# add_if_exists(mpsc-queue)

# add_if_exists(futex)
# add_if_exists(hazard-ptr)
# add_if_exists(lock-free-stack)
# add_if_exists(sync-map)

# add_if_exists(executor)

# if (UNIX AND NOT APPLE)
#   include(tools/cmake/Protobuf.cmake)
#   if (ENABLE_PRIVATE_TESTS)
#     set(FOLLY_DIR "/opt/installed-folly")
#   else()
#     set(FOLLY_DIR "tools/folly")
#   endif()
#   find_package(Boost HINTS ${FOLLY_DIR})
#   find_package(fmt HINTS ${FOLLY_DIR})
#   find_package(folly HINTS ${FOLLY_DIR})
#   find_package(gflags HINTS ${FOLLY_DIR})

#   if (Protobuf_FOUND AND Boost_FOUND AND fmt_FOUND AND folly_FOUND AND gflags_FOUND)
#     message("Networking tasks are enabled")

#     add_if_exists(tools/cactus)

#     add_if_exists(portscan)
#     add_if_exists(portknock)
#     add_if_exists(socks4)
#     add_if_exists(fanout)

#     add_if_exists(socksd)
#     add_if_exists(balancer)

#     add_if_exists(tools/redis)
#   else()
#     message(WARNING "Networking tasks disabled")
#   endif()

# else()
#   message(WARNING "For networking tasks only Linux is supported. "
#     "See instructions for remote toolchain setup."
#   )
# endif()

# add_if_exists(sfinae)
# add_if_exists(compile-eval)
# add_if_exists(constexpr-map)

# add_if_exists(coroutine)
# add_if_exists(concepts)
# add_if_exists(ranges)

# if (LLVM_FOUND)
#   message("clang-fun enabled")
#   add_if_exists(clang-fun)
# else()
#   message(WARNING "clang-fun disabled")
# endif()

# add_if_exists(matrix-2.0)
# add_if_exists(matrix-3.0)
