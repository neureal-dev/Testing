cmake_minimum_required (VERSION 3.10)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable (Concurency "Concurency.cpp" "Concurency.h")

if (MSVC)
	string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	target_compile_options(Concurency PRIVATE /Oy- /W4 /WX )
else()
	target_compile_options(Concurency PRIVATE -fno-omit-frame-pointer -Wall -Wextra -pedantic -Werror -fcoroutines-ts -stdlib=libc++)
endif()

find_package(benchmark CONFIG REQUIRED)
set_target_properties(benchmark::benchmark benchmark::benchmark_main PROPERTIES
	MAP_IMPORTED_CONFIG_MINSIZEREL Release
	MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(Concurency PUBLIC OpenMP::OpenMP_CXX)
endif()

target_link_libraries(Concurency PRIVATE benchmark::benchmark benchmark::benchmark_main)
