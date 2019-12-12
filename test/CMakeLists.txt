cmake_minimum_required(VERSION 3.13)

include(GoogleTest)

find_package(GTest CONFIG REQUIRED)

file(GLOB unit_test_files *.cpp)

add_executable(test_main ${unit_test_files})

if(MSVC)
	string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
	target_compile_options(test_main PRIVATE /Oy- /W4 /WX)
else()
	target_compile_options(test_main PRIVATE -fno-omit-frame-pointer -Wall -Wextra -pedantic -Werror)
endif()

set_target_properties( GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main PROPERTIES
	MAP_IMPORTED_CONFIG_MINSIZEREL Release
	MAP_IMPORTED_CONFIG_RELWITHDEBINFO Release
)

target_link_libraries(test_main PRIVATE GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)

gtest_discover_tests(test_main TEST_PREFIX UT::)

#add_test(allInTest test_main)
