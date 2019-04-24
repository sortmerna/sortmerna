#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 09, 2017
#
# CMake file for sortmerna/src/indexdb
#===============================================================================

cmake_minimum_required(VERSION 3.13)
project(indexdb CXX)

set(CMAKE_CXX_STANDARD 17)

if(WIN32)
	set(IXDB_HDRS
		${DIRENTWIN_HOME}/include/dirent.h
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
		${IXDB_HDRS} ${IXDB_SRCS}
	)
	target_include_directories(indexdb PUBLIC 
		${CMAKE_SOURCE_DIR}/include
		${DIRENTWIN_HOME}/include
		$<TARGET_PROPERTY:cmph,INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
	)
	add_dependencies(indexdb cmph winapi build_version)
else()
	add_executable(indexdb 
		$<TARGET_OBJECTS:cmph>
		$<TARGET_OBJECTS:build_version>
		${IXDB_HDRS} ${IXDB_SRCS}
	)
	target_include_directories(indexdb PUBLIC 
		${CMAKE_SOURCE_DIR}/include
		$<TARGET_PROPERTY:cmph,INCLUDE_DIRECTORIES>
	)
	add_dependencies(indexdb cmph build_version)
	
	if(PORTABLE)
		target_link_options(indexdb PRIVATE -static)
	endif(PORTABLE)
endif()

message("LINK_OPTIONS: ${LINK_OPTIONS}")

