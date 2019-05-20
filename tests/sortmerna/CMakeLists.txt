#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 06, 2018 Wed
#===============================================================================

cmake_minimum_required(VERSION 3.13)
project(tests CXX)

set(CMAKE_CXX_STANDARD 17)

message("tests PROJECT_NAME = ${PROJECT_NAME}")
message("tests PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
message("tests CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message("tests CMAKE_CFG_INTDIR = ${CMAKE_CFG_INTDIR}")

set(TEST_SRCS
	kvdb.cpp
	main.cpp
)

add_executable(tests ${TEST_SRCS})
#target_include_directories(tests
#PRIVATE
#		${CMAKE_SOURCE_DIR}/include
#		#${CMAKE_SOURCE_DIR}/src/winapi/include
#)
target_link_libraries(tests
	build_version
	alp
	smr_objs
	winapi
	$<TARGET_OBJECTS:cmph>
)

#add_executable("test_${test}" ${test}.cpp $<TARGET_OBJECTS:smr_objs>)