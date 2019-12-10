cmake_minimum_required(VERSION 3.13)

project(Testing)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

# Include sub-projects.

add_subdirectory ("BinarySearch")
#add_subdirectory ("Concurency")
add_subdirectory ("LUTBenchmark")
add_subdirectory ("LUTSearch")
#add_subdirectory ("Performance")
add_subdirectory("test")

#add_test(allInSSTest test/test_main)
