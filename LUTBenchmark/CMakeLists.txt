# CMakeList.txt : CMake project for Test, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Oy-")

find_package(benchmark CONFIG REQUIRED)

# Add source to this project's executable.
add_executable (LUTBenchmark "LUTBenchmark.cpp" "LUTBenchmark.h")

target_link_libraries(LUTBenchmark PRIVATE benchmark::benchmark benchmark::benchmark_main)

# TODO: Add tests and install targets if needed.
