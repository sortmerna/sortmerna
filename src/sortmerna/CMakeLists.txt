#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/src/sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.13)
project(sortmerna CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Threads REQUIRED)

find_package(RocksDB 
	REQUIRED
  PATHS ${ROCKSDB_HOME}
  NO_DEFAULT_PATH
)

# fix for RocksDBTargets.cmake, which makes both static and shared libraries always dependent on libz.so
if(NOT WIN32)
	if(ZLIB_LIBRARY_RELEASE)
		set_target_properties(RocksDB::rocksdb PROPERTIES
			INTERFACE_LINK_LIBRARIES "${ZLIB_LIBRARY_RELEASE};-lpthread"
		)
	endif()
endif()

include(FindZLIB)
# prevent CONFIG search mode
find_package(ZLIB MODULE REQUIRED)

include(FindRapidJson)
 # prevent CONFIG search mode
find_package(RapidJson MODULE REQUIRED)

set(SMR_SRCS
	alignment.cpp
	bitvector.cpp
	callbacks.cpp
	cmd.cpp
	gzip.cpp
	index.cpp
	kseq_load.cpp
	kvdb.cpp
	options.cpp
	output.cpp
	paralleltraversal.cpp
	processor.cpp
	read.cpp
	reader.cpp
	readstats.cpp
	references.cpp
	refstats.cpp
	ssw.c
	traverse_bursttrie.cpp
	util.cpp
	writer.cpp
)

# SMR Objects - build a separate library to use with Tests
add_library(smr_objs OBJECT ${SMR_SRCS})
target_link_libraries(smr_objs
	PUBLIC
		ZLIB::ZLIB
		RocksDB::rocksdb
		RapidJSON::RapidJSON
		alp
)
if(WIN32)
	target_include_directories(smr_objs 
		PRIVATE
			${CMAKE_SOURCE_DIR}/include
			$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
			${DIRENTWIN_HOME}/include
			${EXTERNAL_DEPS}/concurrentqueue
	)
else()
	target_include_directories(smr_objs PRIVATE
		${CMAKE_SOURCE_DIR}/include
		${EXTERNAL_DEPS}/concurrentqueue
	)
endif()

# SMR
add_executable(sortmerna main.cpp)
if(WIN32)
	target_include_directories(sortmerna 
		PRIVATE
			${CMAKE_SOURCE_DIR}/include
			${DIRENTWIN_HOME}/include
			${EXTERNAL_DEPS}/concurrentqueue
	)
	target_link_libraries(sortmerna
		build_version
		smr_objs
		winapi
		alp
		ZLIB::ZLIB
		RocksDB::rocksdb
		RapidJSON::RapidJSON
		Rpcrt4.lib
		Cabinet.lib # rocksdb:XPRESS::CreateCompressor,Compress,CloseCompressor,CreateDecompressor,Decompress,CloseDecompressor
	)
	add_dependencies(sortmerna 
		winapi
		alp
		build_version
		smr_objs
	)
else()
	target_include_directories(sortmerna 
		PRIVATE
			${EXTERNAL_DEPS}/concurrentqueue
			${CMAKE_SOURCE_DIR}/include
	)
	target_link_libraries(sortmerna
		build_version
		smr_objs
		alp
		ZLIB::ZLIB
		RocksDB::rocksdb
		RapidJSON::RapidJSON
		${CMAKE_THREAD_LIBS_INIT}
	)
	add_dependencies(sortmerna 
		alp 
		smr_objs 
		build_version
	)
endif()
