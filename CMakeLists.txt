#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 06, 2017
#
# Top level CMake file for sortmerna
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(sortmerna CXX C)

option(WITH_OPENMP "build with OpenMP support. On Nix-systems only i.e. No Windows for now" ON)

#include(ExternalProject)
message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}") # 'GNU' on Ubuntu

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
	if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND "${EXTRA_CXX_FLAGS_RELEASE}" STRGREATER "")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${EXTRA_CXX_FLAGS_RELEASE}")
		set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_C_FLAGS_RELEASE}")
	endif()

	if(WITH_OPENMP AND NOT WIN32)
		find_package(OpenMP)
		if(OPENMP_FOUND)
			set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
			set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${OpenMP_CXX_FLAGS}")
		endif(OPENMP_FOUND)
	#	set(HAVE_OPENMP "${OPENMP_FOUND}")
	endif()
endif()

set(EXTERNAL_DEPS ${CMAKE_SOURCE_DIR}/3rdparty)

# download zlib
if(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/zlib)
	execute_process(COMMAND git clone --depth 1 https://github.com/madler/zlib
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
	)
endif(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/zlib)

# download dirent. Only dirent.h used from this project - no need to build
if(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/dirent)
	execute_process(COMMAND git clone --depth 1 https://github.com/tronkko/dirent
		WORKING_DIRECTORY ${EXTERNAL_DEPS}
)
endif(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/dirent)

# mman for Windows
if(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/mman-win32)
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
	# remove unused directories
	file(REMOVE_RECURSE
		${EXTERNAL_DEPS}/mman-win32/.svn
		${EXTERNAL_DEPS}/mman-win32/branches
		${EXTERNAL_DEPS}/mman-win32/tags
		${EXTERNAL_DEPS}/mman-win32/trunk
		${EXTERNAL_DEPS}/mman-win32/wiki
	)
endif(WIN32 AND NOT EXISTS ${EXTERNAL_DEPS}/mman-win32)

# build zlib
if(WIN32)
	add_subdirectory(${EXTERNAL_DEPS}/zlib)
endif(WIN32)

# build cmph
add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/cmph)

# build alp
add_subdirectory(${CMAKE_SOURCE_DIR}/3rdparty/alp)

# build mman-win32
if(WIN32)
	add_subdirectory(${EXTERNAL_DEPS}/mman-win32)
endif(WIN32)

# build winapi
if(WIN32)
	add_subdirectory(${CMAKE_SOURCE_DIR}/src/winapi)
endif(WIN32)

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