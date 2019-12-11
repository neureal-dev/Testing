cmake_minimum_required (VERSION 3.8)

add_executable (LUTBenchmark "LUTBenchmark.cpp" "LUTBenchmark.h")

if (MSVC)
	#string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	target_compile_options(LUTBenchmark PRIVATE /Oy- /WX)
else()
	target_compile_options(LUTBenchmark PRIVATE -fno-omit-frame-pointer -Wall -Wextra -pedantic -Werror)
endif()

find_package(benchmark CONFIG REQUIRED)

set_target_properties(benchmark::benchmark PROPERTIES
	MAP_IMPORTED_CONFIG_MINSIZEREL Release
	MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
)

target_link_libraries(LUTBenchmark PRIVATE benchmark::benchmark)
#benchmark::benchmark_main)
