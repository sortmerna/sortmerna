#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 09, 2017
#
# CMake file for sortmerna/src/indexdb
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(indexdb CXX)

if(WIN32)
	set(IXDB_HDRS
		${CMAKE_SOURCE_DIR}/3rdparty/dirent/include/dirent.h
	)
endif(WIN32)

set(IXDB_SRCS
	indexdb.cpp
)

if(WIN32)
	add_executable(indexdb 
		$<TARGET_OBJECTS:cmph>
		$<TARGET_OBJECTS:winapi>
		$<TARGET_OBJECTS:mmanwin32>
		${IXDB_HDRS} ${IXDB_SRCS}
	)
else(WIN32)
	add_executable(indexdb 
		$<TARGET_OBJECTS:cmph>
		${IXDB_HDRS} ${IXDB_SRCS}
	)
endif(WIN32)

if(WIN32)
	target_include_directories(indexdb PUBLIC 
		${CMAKE_SOURCE_DIR}/include 
		${CMAKE_SOURCE_DIR}/3rdparty/dirent/include
		$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:mmanwin32,INCLUDE_DIRECTORIES>
	)
else(WIN32)
	target_include_directories(indexdb PUBLIC 
		${CMAKE_SOURCE_DIR}/include 
	)
endif(WIN32)

if(WIN32)
	add_dependencies(indexdb winapi mmanwin32)
else(WIN32)
	add_dependencies(indexdb)
endif(WIN32)
