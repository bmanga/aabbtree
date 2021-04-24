find_package(doctest REQUIRED)

add_executable(test test.cpp)
target_link_libraries(test PRIVATE doctest::doctest abt)