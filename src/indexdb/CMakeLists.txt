#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 09, 2017
#
# CMake file for sortmerna/src/indexdb
#===============================================================================

cmake_minimum_required(VERSION 3.5)
project(indexdb CXX)

if(WIN32)
	set(IXDB_HDRS
		${CMAKE_SOURCE_DIR}/3rdparty/dirent/include/dirent.h
	)
endif()

set(IXDB_SRCS
	indexdb.cpp
)

if(WIN32)
	add_executable(indexdb 
		$<TARGET_OBJECTS:cmph>
		$<TARGET_OBJECTS:winapi>
		$<TARGET_OBJECTS:build_version>
#		$<TARGET_OBJECTS:mmanwin32>
		${IXDB_HDRS} ${IXDB_SRCS}
	)
else()
	add_executable(indexdb 
		$<TARGET_OBJECTS:cmph>
		$<TARGET_OBJECTS:build_version>
		${IXDB_HDRS} ${IXDB_SRCS}
	)
endif()

if(WIN32)
	target_include_directories(indexdb PUBLIC 
		${CMAKE_SOURCE_DIR}/include
		${CMAKE_SOURCE_DIR}/3rdparty/dirent/include
		$<TARGET_PROPERTY:cmph,INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
#		$<TARGET_PROPERTY:mmanwin32,INCLUDE_DIRECTORIES>
	)
else()
	target_include_directories(indexdb PUBLIC 
		${CMAKE_SOURCE_DIR}/include
		$<TARGET_PROPERTY:cmph,INCLUDE_DIRECTORIES>
	)
endif()

if(WIN32)
	add_dependencies(indexdb cmph winapi build_version) #mmanwin32
else(WIN32)
	add_dependencies(indexdb cmph build_version)
endif(WIN32)
