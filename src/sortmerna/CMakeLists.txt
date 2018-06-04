#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/src/sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.5)
project(sortmerna CXX)

set(CMAKE_CXX_STANDARD 14)

#if(MSVC)
#  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /${RUNTIME_LIBRARY}d")
#  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /${RUNTIME_LIBRARY}")
#endif()

if(WIN32)
	set(SMR_HDRS
		${CMAKE_SOURCE_DIR}/include
#		${EXTERNAL_DEPS}/dirent/include
	)
else()
	set(SMR_HDRS
		${CMAKE_SOURCE_DIR}/include
	)
endif()

find_package(Threads REQUIRED)

set(SMR_SRCS
	alignment.cpp
	bitvector.cpp
	callbacks.cpp
	cmd.cpp
	gzip.cpp
	index.cpp
	kseq_load.cpp
	main.cpp
#	minoccur.cpp
#	mmap.cpp
	output.cpp
	paralleltraversal.cpp
	processor.cpp
	read.cpp
	reader.cpp
	readstats.cpp
	references.cpp
	refstats.cpp
	ssw.c
#	ssw_example.c
	traverse_bursttrie.cpp
	writer.cpp
)

#message("${SMR_SRCS}")

if(SET_ROCKSDB)
	add_library(rocksdb STATIC IMPORTED)
	set_target_properties(rocksdb PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${ROCKSDB_INCLUDE_DIR})
	if(WIN32)
		set_target_properties(rocksdb PROPERTIES IMPORTED_LOCATION_DEBUG ${ROCKSDB_LIB_DEBUG}/rocksdb.lib)
		set_target_properties(rocksdb PROPERTIES IMPORTED_LOCATION_RELEASE ${ROCKSDB_LIB_RELEASE}/rocksdb.lib)
	else()
		set_target_properties(rocksdb PROPERTIES IMPORTED_LOCATION_DEBUG ${ROCKSDB_LIB_DEBUG}/librocksdb.a)
		set_target_properties(rocksdb PROPERTIES IMPORTED_LOCATION_RELEASE ${ROCKSDB_LIB_RELEASE}/librocksdb.a)
	endif()
endif()

if(SET_ZLIB)
	add_library(z STATIC IMPORTED)
	set_target_properties(z PROPERTIES IMPORTED_LOCATION_DEBUG ${ZLIB_LIB_DEBUG}/zlibstaticd.lib)
	set_target_properties(z PROPERTIES IMPORTED_LOCATION_RELEASE ${ZLIB_LIB_RELEASE}/zlibstatic.lib)
	set_target_properties(z PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${EXTERNAL_DEPS}/zlib/dist/include)
endif()

if(WIN32)
	add_executable(sortmerna
		$<TARGET_OBJECTS:winapi>
		${SMR_SRCS} ${SMR_HDRS}
	)
else()
	add_executable(sortmerna
		${SMR_SRCS} ${SMR_HDRS}
	)
endif()

#target_compile_definitions(sortmerna PRIVATE
#	_CRT_SECURE_NO_WARNINGS
#)

if(WIN32)
	target_include_directories(sortmerna PRIVATE
		$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
		${EXTERNAL_DEPS}/rapidjson/include
		${EXTERNAL_DEPS}/concurrentqueue
		${EXTERNAL_DEPS}/dirent/include
		${CMAKE_SOURCE_DIR}/include
	)
else()
	target_include_directories(sortmerna PRIVATE
		${EXTERNAL_DEPS}/rapidjson/include
		${EXTERNAL_DEPS}/concurrentqueue
		${CMAKE_SOURCE_DIR}/include
	)
endif()

#find_library(ROCKSDB_LIBS
#  NAMES rocksdb
#  HINTS ${EXTERNAL_DEPS}/rocksdb/build
#  PATH_SUFFIXES Release Debug
#  NO_DEFAULT_PATH
#)
#message("ROCKSDB_LIBS = ${ROCKSDB_LIBS}")

#${ZLIB_LIBRARY_RELEASE}
# 'optimized ${ROCKSDB_LIB_DIR}/rocksdb.lib' only used when building Release. DEBUG gives linking errors because rocksdb.lib not found
if(WIN32)
	target_link_libraries(
		sortmerna
		alp
		z
		rocksdb
		Rpcrt4.lib
		Cabinet.lib # rocksdb:XPRESS::CreateCompressor,Compress,CloseCompressor,CreateDecompressor,Decompress,CloseDecompressor
	)
else()
	target_link_libraries(
		sortmerna
		alp
		z
		rocksdb
		${CMAKE_THREAD_LIBS_INIT}
	)
endif()

if(WIN32)
	add_dependencies(sortmerna winapi alp)
else()
	add_dependencies(sortmerna alp)
endif()