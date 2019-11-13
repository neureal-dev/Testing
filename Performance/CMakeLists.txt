# CMakeList.txt : CMake project for Performance, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.14)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


find_package(benchmark CONFIG REQUIRED)


#link_libraries(debug ${benchmark_DEBUG} optimized ${benchmark_RELEASE})

#SET( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} "/NODEFAULTLIBS:LIBCMT")
#SET( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} "/NODEFAULTLIBS:MSVCRT")
#SET( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} "/NODEFAULTLIBS:MSVCPRTD")

#set(CMAKE_CXX_FLAGS_RELEASE "/MD")
#set(CMAKE_CXX_FLAGS_DEBUG "/MTd")



#If (MSVC)
#Message ( STATUS "Setting MSVC switches")
#string(REPLACE "/MTd" "/MDd" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
#set(CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
#Endif (MSVC)

#set(BENCHMARK_ENABLE_LTO OFF)
#set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
#ADD_DEFINITIONS(-DNODEFAULTLIB:library)
#include_directories(${benchmark_INCLUDE_DIRS})

#include_directories(headers)

#if(MSVC)
  #set(CMAKE_EXE_LINKER_FLAGS
      #"${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB")
#endif()

# Add source to this project's executable.
add_executable (Performance "Performance.cpp" "Performance.h")

target_link_libraries(Performance PRIVATE benchmark::benchmark benchmark::benchmark_main)

set_target_properties(benchmark::benchmark benchmark::benchmark_main PROPERTIES
	MAP_IMPORTED_CONFIG_MINSIZEREL Release
	MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
)

#target_link_libraries(Performance benchmark)
#target_link_libraries(Performance ${benchmark_LIBRARIES})


#if ("${CMAKE_SYSTEM_NAME}" MATCHES "Windows")
 #   target_link_libraries(Performance Shlwapi)
#endif()
# TODO: Add tests and install targets if needed.
