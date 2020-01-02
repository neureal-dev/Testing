cmake_minimum_required(VERSION 3.13)

# Include directories for source code and any used libraries
include_directories(${CMAKE_SOURCE_DIR}/src)

# recursively capture all cpp files under test/src for the test runner
file(GLOB_RECURSE benchmark_files CONFIGURE_DEPENDS src/*.cpp)

add_executable(testing_benchmark main.cpp ${benchmark_files})

if (MSVC)
	target_compile_options(testing_benchmark PRIVATE /Oy- /Zi)
else()
	#flags for clang not clang-cl
	target_compile_options(testing_benchmark PRIVATE -fno-omit-frame-pointer)
endif()

find_package(benchmark CONFIG REQUIRED)

set_target_properties(benchmark::benchmark PROPERTIES
	MAP_IMPORTED_CONFIG_MINSIZEREL Release
	MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
)

target_link_libraries(testing_benchmark PRIVATE benchmark::benchmark)
#benchmark::benchmark_main)
