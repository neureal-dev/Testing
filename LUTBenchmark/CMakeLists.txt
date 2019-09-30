# CMakeList.txt : CMake project for Test, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Oy-")

option (ENABLE_FRAME_POINTER "Disable frame pointer omission" ON)

if (ENABLE_FRAME_POINTER) # disable frame pointer omission
    if (MSVC)
	    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Oy-")
    else()
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
    endif()
endif()

find_package(benchmark CONFIG REQUIRED)

# Add source to this project's executable.
add_executable (LUTBenchmark "LUTBenchmark.cpp" "LUTBenchmark.h")

SET_TARGET_PROPERTIES(LUTBenchmark PROPERTIES LINK_FLAGS "/PROFILE")

target_link_libraries(LUTBenchmark PRIVATE benchmark::benchmark benchmark::benchmark_main)



# TODO: Add tests and install targets if needed.
