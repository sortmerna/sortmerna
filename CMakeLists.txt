#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 06, 2017
#
# Top level CMake file for sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(sortmerna CXX C)

#include(ExternalProject)

set(EXTERNAL_DEPS ${CMAKE_SOURCE_DIR}/3rdparty)

# download zlib
if(NOT EXISTS ${EXTERNAL_DEPS}/zlib/.git)
	execute_process(COMMAND git clone --depth 1 https://github.com/madler/zlib
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
	)
endif(NOT EXISTS ${EXTERNAL_DEPS}/zlib/.git)

# download dirent. Only dirent.h used from this project - no need to build
if(NOT EXISTS ${EXTERNAL_DEPS}/dirent/.git)
	execute_process(COMMAND git clone --depth 1 https://github.com/tronkko/dirent
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
)
endif(NOT EXISTS ${EXTERNAL_DEPS}/dirent/.git)

# mman for Windows
if(NOT EXISTS ${EXTERNAL_DEPS}/mman-win32)
	# download archive. The archive contains directory "mman-win32".
	file(DOWNLOAD
		https://storage.googleapis.com/google-code-archive-source/v2/code.google.com/mman-win32/source-archive.zip
		${EXTERNAL_DEPS}/mman.zip
	)
	# extract mman-win32 archive.
	execute_process(COMMAND ${CMAKE_COMMAND} -E tar -xf mman.zip
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
	)
	# remove mman-win32 archive
	file(REMOVE ${EXTERNAL_DEPS}/mman.zip)
	# move source file. NOTE: modify mman.c to include 'sys/mman.h' instead of 'mman.h'
	file(RENAME ${EXTERNAL_DEPS}/mman-win32/trunk/mman.c ${EXTERNAL_DEPS}/mman-win32/mman.c)
	# move header file.
	file(MAKE_DIRECTORY ${EXTERNAL_DEPS}/mman-win32/sys)
	file(RENAME ${EXTERNAL_DEPS}/mman-win32/trunk/mman.h ${EXTERNAL_DEPS}/mman-win32/sys/mman.h)
endif(NOT EXISTS ${EXTERNAL_DEPS}/mman-win32)

# build zlib
add_subdirectory(${EXTERNAL_DEPS}/zlib)

# build cmph
add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/cmph)

# build mman-win32
add_subdirectory(${EXTERNAL_DEPS}/mman-win32)

# build winapi
add_subdirectory(${CMAKE_SOURCE_DIR}/src/winapi)

add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/alp)

# build indexdb
add_subdirectory(${CMAKE_SOURCE_DIR}/src/indexdb)

# build sortmerna
add_subdirectory(${CMAKE_SOURCE_DIR}/src/sortmerna)


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