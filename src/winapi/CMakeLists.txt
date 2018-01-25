#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 10, 2017
#
# CMake file for sortmerna/src/winapi
#===============================================================================

cmake_minimum_required(VERSION 3.5)
project(winapi CXX)

set(WINAPI_HDRS
	include/unistd.h
	include/sys/time.h
)

#add_subdirectory(${CMAKE_SOURCE_DIR}/winapi)

set(WINAPI_SRCS
	unistd.cpp
	wintime.cpp
)

#message("${WINAPI_SRCS}")

add_library(winapi OBJECT ${WINAPI_SRCS} ${WINAPI_HDRS})

source_group(sys FILES ${CMAKE_SOURCE_DIR}/src/winapi/include/sys/time.h) # Not working: "Header Files\\sys" "include\\sys" "sys"

target_include_directories(winapi PUBLIC ${CMAKE_SOURCE_DIR}/src/winapi/include)