# CMakeList.txt : CMake project for Performance, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(benchmark CONFIG REQUIRED)

# Add source to this project's executable.
add_executable (Concurency "Concurency.cpp" "Concurency.h")

target_link_libraries(Concurency PRIVATE benchmark::benchmark benchmark::benchmark_main)
