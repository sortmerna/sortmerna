#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 06, 2017
#
# Top level CMake file for sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.5)
project(sortmerna CXX C)

set(CMAKE_CXX_STANDARD 14)

set(EXTERNAL_DEPS ${CMAKE_SOURCE_DIR}/3rdparty)

#include(ExternalProject)
message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}") # 'GNU' on Ubuntu
message("CMAKE_CONFIGURATION_TYPES = ${CMAKE_CONFIGURATION_TYPES}")
message("CMAKE_HOST_SYSTEM_NAME = ${CMAKE_HOST_SYSTEM_NAME}")
message("CMAKE_HOST_SYSTEM_VERSION = ${CMAKE_HOST_SYSTEM_VERSION}")
message("CMAKE_HOST_SYSTEM = ${CMAKE_HOST_SYSTEM}")

# append sortmerna custom modules to cmake modules
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules/")

option(SRC_ZLIB "Download Zlib sources" OFF)
option(SRC_ROCKSDB "Download RocksDB sources" OFF)
option(SRC_RAPIDJSON "Download Rapidjson sources" OFF)
option(SET_ROCKSDB "Add RocksDB headers and lib to the build (if RocksDB was built from sources)" OFF)
option(SET_ZLIB "Add Zlib headers and lib to the build (if Zlib was built from sources)" OFF)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND "${EXTRA_CXX_FLAGS_RELEASE}" STRGREATER "")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_CXX_FLAGS_RELEASE}")
		set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_C_FLAGS_RELEASE}")
	endif()
endif()

message("CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message("EXTRA_CXX_FLAGS_RELEASE: ${EXTRA_CXX_FLAGS_RELEASE}")

## download zlib
if(SRC_ZLIB AND NOT EXISTS ${EXTERNAL_DEPS}/zlib)
	execute_process(COMMAND git clone --depth 1 https://github.com/madler/zlib
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
	)
endif()

## Download RocksDB
if(SRC_ROCKSDB AND NOT EXISTS ${EXTERNAL_DEPS}/rocksdb)
	execute_process(COMMAND git clone --depth 1 https://github.com/facebook/rocksdb
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
	)
endif()

## Download Rapidjson
if(SRC_RAPIDJSON AND NOT EXISTS ${EXTERNAL_DEPS}/rapidjson)
	execute_process(COMMAND git clone --depth 1 https://github.com/Tencent/rapidjson
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
	)
endif()

# download dirent. Only dirent.h used from this project - no need to build
if(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/dirent)
	execute_process(COMMAND git clone --depth 1 https://github.com/tronkko/dirent
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
)
endif()

## Download concurrentqueue
if(NOT EXISTS ${EXTERNAL_DEPS}/concurrentqueue)
	execute_process(COMMAND git clone --depth 1 https://github.com/cameron314/concurrentqueue
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
	)
endif()

## mman for Windows TODO: not used any longer. Remove.
#if(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/mman-win32)
	## download archive. The archive contains directory "mman-win32".
#	file(DOWNLOAD
#		https://storage.googleapis.com/google-code-archive-source/v2/code.google.com/mman-win32/source-archive.zip
#		${EXTERNAL_DEPS}/mman.zip
#	)
	## extract mman-win32 archive.
#	execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf mman.zip
#		WORKING_DIRECTORY ${EXTERNAL_DEPS}
#	)
	## remove mman-win32 archive
#	file(REMOVE ${EXTERNAL_DEPS}/mman.zip)
	## move source file. NOTE: modify mman.c to include 'sys/mman.h' instead of 'mman.h'
#	file(RENAME ${EXTERNAL_DEPS}/mman-win32/trunk/mman.c ${EXTERNAL_DEPS}/mman-win32/mman.c)
	## move header file.
#	file(MAKE_DIRECTORY ${EXTERNAL_DEPS}/mman-win32/sys)
#	file(RENAME ${EXTERNAL_DEPS}/mman-win32/trunk/mman.h ${EXTERNAL_DEPS}/mman-win32/sys/mman.h)
	## remove unused directories
#	file(REMOVE_RECURSE
#		${EXTERNAL_DEPS}/mman-win32/.svn
#		${EXTERNAL_DEPS}/mman-win32/branches
#		${EXTERNAL_DEPS}/mman-win32/tags
#		${EXTERNAL_DEPS}/mman-win32/trunk
#		${EXTERNAL_DEPS}/mman-win32/wiki
#	)
#endif(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/mman-win32)

## build cmph
add_subdirectory(${EXTERNAL_DEPS}/cmph)

## build alp
add_subdirectory(${EXTERNAL_DEPS}/alp)

## build winapi
if(WIN32)
	add_subdirectory(${CMAKE_SOURCE_DIR}/src/winapi)
endif()

## build indexdb
add_subdirectory(${CMAKE_SOURCE_DIR}/src/indexdb)

## build sortmerna
add_subdirectory(${CMAKE_SOURCE_DIR}/src/sortmerna)

## cache entries for 3rd party libraries
if(SET_ROCKSDB)
	## cache entries for RocksDB
	set(ROCKSDB_INCLUDE_DIR ${EXTERNAL_DEPS}/rocksdb/include CACHE PATH "Rocksdb Include directory")
	set(ROCKSDB_LIB_DEBUG ${EXTERNAL_DEPS}/rocksdb/build/Debug CACHE PATH "Rocksdb Debug Lib directory")
	set(ROCKSDB_LIB_RELEASE ${EXTERNAL_DEPS}/rocksdb/build/Release CACHE PATH "Rocksdb Release Lib directory")
	message("ROCKSDB_INCLUDE_DIR = ${ROCKSDB_INCLUDE_DIR}")
	message("CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
endif()

if(SET_ZLIB)
	## zlib
	set(ZLIB_LIB_DEBUG ${EXTERNAL_DEPS}/zlib/build/Debug CACHE PATH "Zlib library directory")
	set(ZLIB_LIB_RELEASE ${EXTERNAL_DEPS}/zlib/build/Release CACHE PATH "Zlib library directory")
endif()

#find_package(zlib REQUIRED)

#option(WITH_MD_LIBRARY "build with MD" ON)
#if(WIN32 AND MSVC)
#  if(WITH_MD_LIBRARY)
#    set(RUNTIME_LIBRARY "MD")
#  else()
#    set(RUNTIME_LIBRARY "MT")
#  endif()
#endif()

#ExternalProject_Add(zlib
#	DOWNLOAD_COMMAND C:/a03_libs/git-2.12.2-x64/bin/git clone --depth 1 https://github.com/madler/zlib
#	UPDATE_COMMAND C:/a03_libs/git-2.12.2-x64/bin/git pull
#	DOWNLOAD_DIR ${EXTERNAL_DEPS_ROOT}  # seems DOWNLOAD_DIR only works with DOWNLOAD_COMMAND. Ignored with GIT_REPOSITORY
#    GIT_REPOSITORY https://github.com/madler/zlib
#	GIT_SHALLOW 1
#    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_DEPS_ROOT}
#)

#include_directories(${EXTERNAL_DEPS_ROOT}/include)
#link_directories(${EXTERNAL_DEPS_ROOT}/lib)

#add_executable(FirstTest main.cpp)
#add_dependencies(FirstTest googletest)
#target_link_libraries(FirstTest gtest gtest_main pthread)