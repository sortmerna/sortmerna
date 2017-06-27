#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/src/sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(sortmerna CXX)

if(WIN32)
	set(SMR_HDRS
		${CMAKE_SOURCE_DIR}/include
		${CMAKE_SOURCE_DIR}/3rdparty/dirent/include/dirent.h
	)
else(WIN32)
	set(SMR_HDRS
		${CMAKE_SOURCE_DIR}/include
	)
endif(WIN32)

set(SMR_SRCS
	alignment.cpp
	bitvector.cpp
	kseq_load.cpp
	load_index.cpp
	main.cpp
#	minoccur.cpp
	mmap.cpp
	outputformats.cpp
	paralleltraversal.cpp
	ssw.c
#	ssw_example.c
	traverse_bursttrie.cpp
)

#message("${SMR_SRCS}")

if(WIN32)
	add_executable(sortmerna
		$<TARGET_OBJECTS:winapi>
		$<TARGET_OBJECTS:mmanwin32> 
		${SMR_SRCS} ${SMR_HDRS}
	)
else(WIN32)
	add_executable(sortmerna
		${SMR_SRCS} ${SMR_HDRS}
	)
endif(WIN32)

#target_compile_definitions(sortmerna PRIVATE
#	_CRT_SECURE_NO_WARNINGS
#)

if(WIN32)
	target_include_directories(sortmerna PRIVATE
		$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:mmanwin32,INCLUDE_DIRECTORIES>
		$<TARGET_PROPERTY:zlib,INCLUDE_DIRECTORIES>
		${CMAKE_SOURCE_DIR}/3rdparty/dirent/include
		${CMAKE_SOURCE_DIR}/include
	)
else(WIN32)
	target_include_directories(sortmerna PRIVATE
		${CMAKE_SOURCE_DIR}/include
	)
endif(WIN32)

if(WIN32)
	target_link_libraries(sortmerna
		alp
		zlibstatic
	)
else(WIN32)
	target_link_libraries(sortmerna
		alp
	)
endif(WIN32)

if(WIN32)
	add_dependencies(sortmerna winapi mmanwin32 zlib alp)
else(WIN32)
	add_dependencies(sortmerna alp)
endif(WIN32)