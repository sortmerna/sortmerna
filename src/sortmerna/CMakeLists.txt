#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/src/sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.13)
project(sortmerna CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(Threads REQUIRED)

find_package(RocksDB REQUIRED
  PATHS ${ROCKSDB_HOME}
  NO_DEFAULT_PATH
)

get_property(trans_deps TARGET RocksDB::rocksdb PROPERTY INTERFACE_LINK_LIBRARIES)
message("rocksdb transitive link dependencies: ${trans_deps}")

# fix RocksDBTargets.cmake, which makes both static and shared libraries always dependent on libz.so
if(NOT WIN32)
	if(ZLIB_LIBRARY_RELEASE)
		# modify the dependent libraries of the RocksDB::rocksdb target
		set_target_properties(RocksDB::rocksdb PROPERTIES
			INTERFACE_LINK_LIBRARIES "${ZLIB_LIBRARY_RELEASE};Threads::Threads"
		)
		get_property(trans_deps TARGET RocksDB::rocksdb PROPERTY INTERFACE_LINK_LIBRARIES)
		message("rocksdb transitive link dependencies AFTER modification: ${trans_deps}")
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
	read_control.cpp
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
		alp
		RocksDB::rocksdb
		RapidJSON::RapidJSON
		ZLIB::ZLIB
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

get_property(trans_deps TARGET smr_objs PROPERTY INTERFACE_LINK_LIBRARIES)
message("SMR Objects transitive link dependencies: ${trans_deps}")

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
else()
	target_include_directories(sortmerna 
		PRIVATE
			${EXTERNAL_DEPS}/concurrentqueue
			${CMAKE_SOURCE_DIR}/include
	)
	#get_property(rocksdb_trans_deps TARGET RocksDB::rocksdb PROPERTY INTERFACE_LINK_LIBRARIES)
	#message("rocksdb_trans_deps: ${rocksdb_trans_deps}")
	target_link_libraries(sortmerna
		build_version
		smr_objs
		alp # this is a transitive dependency of smr_objs but the linking fails without it. Why?
		# the following are all transitive dependencies of smr_objs i.e. no need to link: 
		# RocksDB::rocksdb RapidJSON::RapidJSON ZLIB::ZLIB Threads::Threads (rockdb deps)
	)
	
	if(PORTABLE)
		target_link_options(sortmerna PRIVATE -static)
	endif(PORTABLE)
endif()
