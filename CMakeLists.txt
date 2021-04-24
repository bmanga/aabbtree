project(aabbtree)

add_library(abt INTERFACE include/abt/aabb_tree.hpp)
target_include_directories(abt INTERFACE include/)

add_executable(demo demos/hard_disc.cc)
target_link_libraries(demo abt)
