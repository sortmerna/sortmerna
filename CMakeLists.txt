#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 06, 2017
#
# Top level CMake file for sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.21)
project(sortmerna CXX C)

set(CMAKE_CXX_STANDARD 17)

find_package(Git)

string(TIMESTAMP GIT_DATE_TIME "%Y/%m/%d %H:%M:%S" UTC)

#
# Prepare Version information
#
if(GIT_FOUND AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/.git")
  if(WIN32)
    execute_process(COMMAND $ENV{COMSPEC} /C ${GIT_EXECUTABLE} -C ${CMAKE_CURRENT_SOURCE_DIR} rev-parse HEAD OUTPUT_VARIABLE GIT_SHA)
  else()
    execute_process(COMMAND ${GIT_EXECUTABLE} -C ${CMAKE_CURRENT_SOURCE_DIR} rev-parse HEAD OUTPUT_VARIABLE GIT_SHA)
  endif()
else()
  set(GIT_SHA 0)
endif()

string(REGEX REPLACE "[^0-9a-f]+" "" GIT_SHA "${GIT_SHA}")

# Read sortmerna version from version.h file.
file(READ include/version.h version_header_file)
string(REGEX MATCH "#define SORTMERNA_MAJOR ([0-9]+)" _ ${version_header_file})
set(SORTMERNA_VERSION_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "#define SORTMERNA_MINOR ([0-9]+)" _ ${version_header_file})
set(SORTMERNA_VERSION_MINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "#define SORTMERNA_PATCH ([0-9]+)" _ ${version_header_file})
set(SORTMERNA_VERSION_PATCH ${CMAKE_MATCH_1})
set(SORTMERNA_VERSION ${SORTMERNA_VERSION_MAJOR}.${SORTMERNA_VERSION_MINOR}.${SORTMERNA_VERSION_PATCH})


#include(ExternalProject)
message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}") # 'GNU' on Ubuntu
message("CMAKE_CONFIGURATION_TYPES = ${CMAKE_CONFIGURATION_TYPES}")
message("CMAKE_HOST_SYSTEM_NAME = ${CMAKE_HOST_SYSTEM_NAME}")
message("CMAKE_HOST_SYSTEM_VERSION = ${CMAKE_HOST_SYSTEM_VERSION}")
message("CMAKE_HOST_SYSTEM = ${CMAKE_HOST_SYSTEM}")
message("PROJECT_BINARY_DIR = ${PROJECT_BINARY_DIR}")
message("CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")
message("CMAKE_CURRENT_BINARY_DIR = ${CMAKE_CURRENT_BINARY_DIR}")
message("RUNTIME_OUTPUT_DIRECTORY_RELEASE = ${RUNTIME_OUTPUT_DIRECTORY_RELEASE}")

# append sortmerna custom modules to cmake modules
message("CMAKE_CURRENT_LIST_DIR: ${CMAKE_CURRENT_LIST_DIR}")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules/")

set(EXTERNAL_DEPS ${CMAKE_SOURCE_DIR}/3rdparty CACHE PATH "Installation directory for 3rd party dependencies")
set(ZLIB_URL https://github.com/madler/zlib.git CACHE STRING "ZLib GitHub URL")
#set(ZLIB_ROOT ${EXTERNAL_DEPS}/zlib CACHE PATH "ZLIB installation directory")
set(ROCKSDB_URL https://github.com/facebook/rocksdb.git CACHE STRING "RocksDB GitHub URL")
set(ROCKSDB_SRC ${EXTERNAL_DEPS}/rocksdb CACHE PATH "RocksDB source root directory")
set(ROCKSDB_DIST ${EXTERNAL_DEPS}/rocksdb/dist CACHE PATH "RocksDB installation directory")
#set(RAPIDJSON_URL https://github.com/Tencent/rapidjson CACHE STRING "RapidJSON GitHub URL")
#set(RAPIDJSON_HOME ${EXTERNAL_DEPS}/rocksdb CACHE PATH "RapidJSON installation directory")
set(DIRENTWIN_URL https://github.com/tronkko/dirent CACHE STRING "Dirent-Windows interface GitHub URL")
set(DIRENTWIN_HOME ${EXTERNAL_DEPS}/dirent CACHE PATH "Dirent-Windows interface installation directory")
set(CONCURRENTQUEUE_URL https://github.com/cameron314/concurrentqueue CACHE STRING "Concurrent queue library GitHub URL")
set(CONCURRENTQUEUE_HOME ${EXTERNAL_DEPS}/concurrentqueue CACHE PATH "Concurrent queue library installation directory")
#set(gflags_DIR  ${EXTERNAL_DEPS}/gflags/dist/lib/cmake/gflags CACHE PATH "GFlags library configuration")

option(ROCKSDB_USE_STATIC_LIBS "Use static RocksDB library" ON)
option(ZLIB_STATIC "Use static ZLib library" ON)

message("CMAKE_CXX_FLAGS_DEBUG: ${CMAKE_CXX_FLAGS_DEBUG}")
message("CMAKE_C_FLAGS_DEBUG: ${CMAKE_C_FLAGS_DEBUG}")

option(WITH_TESTS "Select whether to build tests" OFF)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND "${EXTRA_CXX_FLAGS_RELEASE}" STRGREATER "")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_CXX_FLAGS_RELEASE}")
		set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_C_FLAGS_RELEASE}")
	endif()
endif()

message("CMAKE_CXX_FLAGS_RELEASE: ${CMAKE_CXX_FLAGS_RELEASE}")
message("EXTRA_CXX_FLAGS_RELEASE: ${EXTRA_CXX_FLAGS_RELEASE}")

## download zlib
#if(SRC_ZLIB AND NOT EXISTS ${EXTERNAL_DEPS}/zlib)
#	execute_process(COMMAND git clone --depth 1 https://github.com/madler/zlib
#		WORKING_DIRECTORY ${EXTERNAL_DEPS}
#	)
#endif()

## Download RocksDB
#if(SRC_ROCKSDB AND NOT EXISTS ${EXTERNAL_DEPS}/rocksdb)
#	execute_process(COMMAND git clone --depth 1 https://github.com/facebook/rocksdb
#		WORKING_DIRECTORY ${EXTERNAL_DEPS}
#	)
#endif()

#
# modify RocksDB thirdparty.inc - no need, it is modified before building RocksDB
#
#option(ROCKSDB_3RDPARTY_INC "Modify rocksdb/thirdparty.inc" OFF)
#if (ROCKSDB_3RDPARTY_INC)
#	include(smr_utils)
#	rocksdb_edit_3rdparty_inc()
#endif(ROCKSDB_3RDPARTY_INC)


## Download Rapidjson
#(SRC_RAPIDJSON AND NOT EXISTS ${EXTERNAL_DEPS}/rapidjson)
#	execute_process(COMMAND git clone --depth 1 https://github.com/Tencent/rapidjson
#		WORKING_DIRECTORY ${EXTERNAL_DEPS}
#	)
#endif()

# download dirent. Only dirent.h used from this project - no need to build
if(WIN32 AND NOT EXISTS ${DIRENTWIN_HOME})
	execute_process(COMMAND git clone --depth 1 ${DIRENTWIN_URL}
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
)
endif()

## Download concurrentqueue
#if(NOT EXISTS ${EXTERNAL_DEPS}/concurrentqueue)
#	execute_process(COMMAND git clone --depth 1 https://github.com/cameron314/concurrentqueue
#		WORKING_DIRECTORY ${EXTERNAL_DEPS}
#	)
#endif()

# build build_version
set(BUILD_VERSION_CPP ${CMAKE_BINARY_DIR}/build_version.cpp)
configure_file(src/util/build_version.cpp.in ${BUILD_VERSION_CPP} @ONLY)
add_library(build_version OBJECT ${BUILD_VERSION_CPP})
target_include_directories(build_version PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

## build cmph
add_subdirectory(${EXTERNAL_DEPS}/cmph)

## build alp
add_subdirectory(${EXTERNAL_DEPS}/alp)

## build winapi
if(WIN32)
	add_subdirectory(${CMAKE_SOURCE_DIR}/src/winapi)
endif()

## build indexdb
#add_subdirectory(${CMAKE_SOURCE_DIR}/src/indexdb)

## build sortmerna
add_subdirectory(${CMAKE_SOURCE_DIR}/src/sortmerna)

if (WITH_TESTS)
  add_subdirectory (tests/sortmerna)
endif ()

# Installation and packaging
#set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/dist CACHE PATH "Install path prefix, prepended onto install directories.")

#include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(package_config_destination ${CMAKE_INSTALL_PREFIX}/cmake/sortmerna)
message("package_config_destination: ${package_config_destination}")

configure_package_config_file(
	${CMAKE_CURRENT_LIST_DIR}/cmake/sortmernaConfig.cmake.in sortmernaConfig.cmake
	INSTALL_DESTINATION ${package_config_destination}
)

write_basic_package_version_file(
	sortmernaConfigVersion.cmake
	VERSION ${SORTMERNA_VERSION}
	COMPATIBILITY SameMajorVersion
)  

install(TARGETS sortmerna
	EXPORT sortmernaTargets
	RUNTIME 
	DESTINATION bin
	COMPONENT runtime
)

install(FILES
		${CMAKE_CURRENT_BINARY_DIR}/sortmernaConfig.cmake
		${CMAKE_CURRENT_BINARY_DIR}/sortmernaConfigVersion.cmake
	DESTINATION ${package_config_destination}
	COMPONENT config
)

# generate PACKAGE and PACKAGE_SOURCE targets, and CPack configuration files
# to build distribution archive/installer
set(CPACK_PACKAGE_VENDOR "Clarity Genomics")
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_INSTALL_PREFIX})
set(CPACK_PACKAGE_VERSION_MAJOR ${SORTMERNA_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${SORTMERNA_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${SORTMERNA_VERSION_PATCH})
if(UNIX AND NOT APPLE) # CMAKE_HOST_SYSTEM_NAME = Linux
	set(CPACK_GENERATOR "DEB" "STGZ" "TGZ") # "TGZ" - use for Conda build to pick-up
	set(CPACK_DEBIAN_PACKAGE_DEBUG ON) # output Deb generator trace debug information
	#set(CPACK_DEB_COMPONENT_INSTALL ON) # Enable component packaging for CPackDEB
	if(NOT PORTABLE)
		set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
	endif()
	set(CPACK_DEBIAN_PACKAGE_MAINTAINER "biocodz")
endif()
include(CPack)