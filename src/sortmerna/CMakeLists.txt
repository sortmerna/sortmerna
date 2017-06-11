#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/src/sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(sortmerna CXX)

set(SMR_HDRS
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/3rdparty/dirent/include/dirent.h
)

set(SMR_SRCS
	${CMAKE_SOURCE_DIR}/src/alignment.cpp
	${CMAKE_SOURCE_DIR}/src/bitvector.cpp
#	indexdb.cpp
	${CMAKE_SOURCE_DIR}/src/kseq_load.cpp
	${CMAKE_SOURCE_DIR}/src/load_index.cpp
	${CMAKE_SOURCE_DIR}/src/main.cpp
#	minoccur.cpp
	${CMAKE_SOURCE_DIR}/src/mmap.cpp
	${CMAKE_SOURCE_DIR}/src/outputformats.cpp
	${CMAKE_SOURCE_DIR}/src/paralleltraversal.cpp
	${CMAKE_SOURCE_DIR}/src/ssw.c
#	ssw_example.c
	${CMAKE_SOURCE_DIR}/src/traverse_bursttrie.cpp
)

#message("${SMR_SRCS}")

add_executable(sortmerna
	$<TARGET_OBJECTS:winapi>
	$<TARGET_OBJECTS:mmanwin32> 
	${SMR_SRCS} ${SMR_HDRS}
)

#target_compile_definitions(sortmerna PRIVATE
#	_CRT_SECURE_NO_WARNINGS
#)

target_include_directories(sortmerna PRIVATE
	$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
	$<TARGET_PROPERTY:mmanwin32,INCLUDE_DIRECTORIES>
	$<TARGET_PROPERTY:zlib,INCLUDE_DIRECTORIES>
	${CMAKE_SOURCE_DIR}/3rdparty/dirent/include
	${CMAKE_SOURCE_DIR}/include
)

#TODO: link zlibstaticd.lib properly
target_link_libraries(sortmerna
	alp
	${CMAKE_SOURCE_DIR}/build/3rdparty/zlib/Debug/zlibstaticd.lib
)