project(aabbtree)
cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 20)
cmake_policy(SET CMP0069 NEW)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)

find_package(Vc REQUIRED)
add_library(abt INTERFACE)
add_compile_options("-Ofast -flto")
target_include_directories(abt INTERFACE include/)
target_link_libraries(abt INTERFACE Vc::Vc)

add_executable(demo demos/hard_disc.cc)
target_link_libraries(demo abt)

set_property(TARGET demo PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)

option(BUILD_TESTS "Build the unit tests" OFF)

if (BUILD_TESTS)
    message(STATUS "Building tests...")
    add_subdirectory(test/)
endif()
