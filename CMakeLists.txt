cmake_minimum_required(VERSION 3.13)

project(Testing)

option(PEDANTIC "Build with high compiler warning level" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (PEDANTIC)
	if (MSVC)
		string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	else()
		compile_options(-Wall -Wextra -pedantic -Werror)
	endif()
endif()

enable_testing()

# Include sub-projects.

add_subdirectory("benchmark")
#add_subdirectory ("BinarySearch")
#add_subdirectory ("Concurency")
add_subdirectory ("LUTBenchmark")
#add_subdirectory ("LUTSearch")
#add_subdirectory ("Performance")
#add_subdirectory("test")

#add_test(allInSSTest test/test_main)
