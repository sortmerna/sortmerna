#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/src/sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.12)
project(sortmerna CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Threads REQUIRED)

find_package(RocksDB 
	REQUIRED
#  5.18.0
#  NAMES rocksdb
  PATHS ${ROCKSDB_HOME}
#  PATH_SUFFIXES Release Debug
  NO_DEFAULT_PATH
)
#message("ROCKSDB_LIBS = ${ROCKSDB_LIBS}")

include(FindZLIB)
find_package(ZLIB 
	REQUIRED
#	PATHS ${ZLIB_HOME} # BUG: if specified gives error 'missing ZLIB_DIR', otherwise works fine.
#	NO_DEFAULT_PATH
)

find_package(RapidJson)

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
)
if(WIN32)
	target_include_directories(smr_objs 
		PRIVATE
			${CMAKE_SOURCE_DIR}/include
			$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
			${EXTERNAL_DEPS}/dirent/include
			${EXTERNAL_DEPS}/concurrentqueue
	)
else()
	target_include_directories(smr_objs PRIVATE
		${CMAKE_SOURCE_DIR}/include
		${EXTERNAL_DEPS}/concurrentqueue
	)
endif()

# SMR
if(WIN32)
	add_executable(sortmerna
		main.cpp
		#$<TARGET_OBJECTS:winapi> # No more: CMake 3.12 allows linking OJBECT libraries
		#$<TARGET_OBJECTS:smr_objs>
		#$<TARGET_OBJECTS:build_version>
	)
	target_include_directories(sortmerna 
		PRIVATE
			${CMAKE_SOURCE_DIR}/include
			${EXTERNAL_DEPS}/concurrentqueue
			${EXTERNAL_DEPS}/dirent/include
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
		winapi alp 
		smr_objs 
		build_version
	)
else()
	add_executable(sortmerna
		#$<TARGET_OBJECTS:smr_objs>
		$<TARGET_OBJECTS:build_version>
		main.cpp
	)
	target_include_directories(sortmerna 
		PRIVATE
			${EXTERNAL_DEPS}/concurrentqueue
			${CMAKE_SOURCE_DIR}/include
	)
	target_link_libraries(sortmerna
		alp
		#z
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

    #include(GetPrerequisites)
    #get_prerequisites(sortmerna SMR_DEPS 0 1 ) # executable has to be specified (not sortmerna target)
	#get_prerequisites(sortmerna SMR_DEPS)
    #message(STATUS "SMR_DEPS: ${SMR_DEPS}")

	#get_target_property(SMR_DEPS sortmerna INTERFACE_LINK_LIBRARIES)
	#message(STATUS "INTERFACE_LINK_LIBRARIES: ${SMR_DEPS}")
	#get_target_property(SMR_DEPS sortmerna IMPORTED_LINK_DEPENDENT_LIBRARIES)
	#message(STATUS "IMPORTED_LINK_DEPENDENT_LIBRARIES: ${SMR_DEPS}")
	#get_target_property(SMR_DEPS sortmerna IMPORTED_LINK_INTERFACE_LIBRARIES)
	#message(STATUS "IMPORTED_LINK_INTERFACE_LIBRARIES: ${SMR_DEPS}")
	#get_target_property(SMR_DEPS sortmerna LINK_DEPENDS)
	#message(STATUS "IMPORTED_LINK_INTERFACE_LIBRARIES: ${SMR_DEPS}")
endif()
