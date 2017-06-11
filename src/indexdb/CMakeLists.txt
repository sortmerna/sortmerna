#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 09, 2017
#
# CMake file for sortmerna/src/indexdb
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(indexdb CXX)

set(IXDB_HDRS
	${CMAKE_SOURCE_DIR}/3rdparty/dirent/include/dirent.h
)

set(IXDB_SRCS
	${CMAKE_SOURCE_DIR}/src/indexdb.cpp
#	minoccur.cpp
#	ssw_example.c
)

add_executable(indexdb 
	$<TARGET_OBJECTS:cmph>
	$<TARGET_OBJECTS:winapi>
	$<TARGET_OBJECTS:mmanwin32>
	${IXDB_HDRS} ${IXDB_SRCS}
)

target_include_directories(indexdb PUBLIC 
	${CMAKE_SOURCE_DIR}/include 
	${CMAKE_SOURCE_DIR}/3rdparty/dirent/include
	$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
	$<TARGET_PROPERTY:mmanwin32,INCLUDE_DIRECTORIES>
)
