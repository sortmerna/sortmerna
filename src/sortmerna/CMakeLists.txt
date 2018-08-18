#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/src/sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.5)
project(sortmerna CXX)

set(CMAKE_CXX_STANDARD 14)

find_package(Threads REQUIRED)

#find_library(ROCKSDB_LIBS
#  NAMES rocksdb
#  HINTS ${EXTERNAL_DEPS}/rocksdb/build
#  PATH_SUFFIXES Release Debug
#  NO_DEFAULT_PATH
#)
#message("ROCKSDB_LIBS = ${ROCKSDB_LIBS}")

set(SMR_SRCS
	alignment.cpp
	bitvector.cpp
	callbacks.cpp
	cmd.cpp
	gzip.cpp
	index.cpp
	kseq_load.cpp
	kvdb.cpp
#	main.cpp
#	minoccur.cpp
#	mmap.cpp
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
#	ssw_example.c
	traverse_bursttrie.cpp
	util.cpp
	writer.cpp
)

# SMR Objects - build a separate library to use with Tests
add_library(smr_objs OBJECT ${SMR_SRCS})
if(WIN32)
	target_include_directories(smr_objs PRIVATE
		$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
		${EXTERNAL_DEPS}/rapidjson/include
		${EXTERNAL_DEPS}/concurrentqueue
		${EXTERNAL_DEPS}/dirent/include
		${CMAKE_SOURCE_DIR}/include
		${EXTERNAL_DEPS}/zlib/dist/include
		${ROCKSDB_INCLUDE_DIR}
	)
else()
	target_include_directories(smr_objs PRIVATE
		${EXTERNAL_DEPS}/rapidjson/include
		${EXTERNAL_DEPS}/concurrentqueue
		${CMAKE_SOURCE_DIR}/include
		${EXTERNAL_DEPS}/zlib/dist/include
		${ROCKSDB_INCLUDE_DIR}
	)
endif()

# SMR
if(WIN32)
	add_executable(sortmerna
		$<TARGET_OBJECTS:winapi>
		$<TARGET_OBJECTS:smr_objs>
		main.cpp
	)
	target_include_directories(sortmerna PRIVATE
		$<TARGET_PROPERTY:winapi,INCLUDE_DIRECTORIES>
		${EXTERNAL_DEPS}/rapidjson/include
		${EXTERNAL_DEPS}/concurrentqueue
		${EXTERNAL_DEPS}/dirent/include
		${CMAKE_SOURCE_DIR}/include
	)
	target_link_libraries( sortmerna
		alp
		z
		rocksdb
		Rpcrt4.lib
		Cabinet.lib # rocksdb:XPRESS::CreateCompressor,Compress,CloseCompressor,CreateDecompressor,Decompress,CloseDecompressor
	)
	add_dependencies(sortmerna winapi alp smr_objs)
else()
	add_executable(sortmerna
		$<TARGET_OBJECTS:smr_objs>
		main.cpp
	)
	target_include_directories(sortmerna PRIVATE
		${EXTERNAL_DEPS}/rapidjson/include
		${EXTERNAL_DEPS}/concurrentqueue
		${CMAKE_SOURCE_DIR}/include
	)
	target_link_libraries( sortmerna
		alp
		z
		rocksdb
		${CMAKE_THREAD_LIBS_INIT}
	)
	add_dependencies(sortmerna alp smr_objs)

    #include(GetPrerequisites)
    #get_prerequisites(sortmerna SMR_DEPS 0 1 ) # executable has to be specified (not sortmerna target)
    #message(STATUS "SMR_DEPS: ${SMR_DEPS}")

	get_prerequisits(sortmerna SMR_DEPS)
	message(STATUS ${SMR_DEPS})

	get_target_property(SMR_DEPS sortmerna INTERFACE_LINK_LIBRARIES)
	message(STATUS "INTERFACE_LINK_LIBRARIES: ${SMR_DEPS}")
	get_target_property(SMR_DEPS sortmerna IMPORTED_LINK_DEPENDENT_LIBRARIES)
	message(STATUS "IMPORTED_LINK_DEPENDENT_LIBRARIES: ${SMR_DEPS}")
	get_target_property(SMR_DEPS sortmerna IMPORTED_LINK_INTERFACE_LIBRARIES)
	message(STATUS "IMPORTED_LINK_INTERFACE_LIBRARIES: ${SMR_DEPS}")
	get_target_property(SMR_DEPS sortmerna LINK_DEPENDS)
	message(STATUS "IMPORTED_LINK_INTERFACE_LIBRARIES: ${SMR_DEPS}")
	
	
endif()
