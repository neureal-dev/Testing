# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.8)

project ("Testing")

# Include sub-projects.
add_subdirectory ("ASIO")
add_subdirectory ("Performance")
add_subdirectory ("BinarySearch")
add_subdirectory ("LUTBenchmark")
add_subdirectory ("LUTSearch")
