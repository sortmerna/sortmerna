#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 06, 2018 Wed
#===============================================================================

cmake_minimum_required(VERSION 3.5)
project(tests CXX)

set(CMAKE_CXX_STANDARD 14)

message("tests PROJECT_NAME = ${PROJECT_NAME}")
message("tests PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
message("tests CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
message("tests CMAKE_CFG_INTDIR = ${CMAKE_CFG_INTDIR}")

set(tests
	kvdb
)

include_directories("${CMAKE_SOURCE_DIR}/include" ${CMAKE_SOURCE_DIR}/src/winapi/include)

foreach(test ${tests})
	add_executable("test_${test}" ${test}.cpp $<TARGET_OBJECTS:smr_objs>)
	if(WIN32)
		target_link_libraries( "test_${test}"
			alp
			z
			rocksdb
			Rpcrt4.lib
			Cabinet.lib
		)
	else()
		target_link_libraries( "test_${test}"
			alp
			z
			rocksdb
		)
	endif()
endforeach()

add_dependencies(test_kvdb sortmerna)