#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/3rdparty/mman-win32
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(mmanwin32 C)

set(MMAN_HDRS sys/mman.h)

set(MMAN_SRCS mman.c)

add_library(mmanwin32 OBJECT ${MMAN_SRCS} ${MMAN_HDRS})
target_include_directories(mmanwin32 PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/mman-win32)