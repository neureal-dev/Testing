﻿# CMakeList.txt : CMake project for BinarySearch, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.14)

# Add source to this project's executable.
add_executable (BinarySearch "BinarySearch.cpp" "BinarySearch.h")
target_compile_options(BinarySearch PRIVATE /await)

# TODO: Add tests and install targets if needed.
