cmake_minimum_required(VERSION 3.0)

project("Chan3D")

enable_testing()

# If you use ./googletest/* locally
# add_subdirectory(googletest)

add_executable(Chan3D chan_tests.cpp)

# If you use ./googletest/* locally
# target_link_libraries(Chan3D gtest_main gmock_main)
target_link_libraries(Chan3D gtest gmock pthread)
