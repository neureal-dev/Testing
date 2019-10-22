cmake_minimum_required (VERSION 3.8)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable (LUTBenchmark "LUTBenchmark.cpp" "LUTBenchmark.h")

if (MSVC)
	target_compile_options(LUTBenchmark PRIVATE /W4 /WX)
else()
	target_compile_options(LUTBenchmark PRIVATE -Wall -Wextra -pedantic -Werror)
endif()

find_package(benchmark CONFIG REQUIRED)

set_target_properties(benchmark::benchmark benchmark::benchmark_main PROPERTIES
	MAP_IMPORTED_CONFIG_MINSIZEREL Release
	MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
)

target_link_libraries(LUTBenchmark PRIVATE benchmark::benchmark benchmark::benchmark_main)
