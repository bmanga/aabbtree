project(aabbtree)

add_library(aabbtree src/AABB.cc "include/aabb/AABB.h")
target_include_directories(aabbtree PUBLIC include/)

add_executable(demo demos/hard_disc.cc)
target_link_libraries(demo aabbtree)
