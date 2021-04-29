project(aabbtree)
cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
add_library(abt INTERFACE include/abt/aabb_tree.hpp)
target_include_directories(abt INTERFACE include/)

add_executable(demo demos/hard_disc.cc)
target_link_libraries(demo abt)

option(BUILD_TESTS "Build the unit tests" OFF)

if (BUILD_TESTS)
    message(STATUS "Building tests...")
    add_subdirectory(test/)
endif()
