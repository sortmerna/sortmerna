#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/cmph
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(cmph CXX)

#AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/cmph CMPH_SRCS)
#include_directories(${CMAKE_SOURCE_DIR}/cmph)              # not recommended practice. 
#file(GLOB CMPH_HDRS ${CMAKE_SOURCE_DIR}/cmph "*.h" )       # not recommended practice. List files explicitly

# find C:/a03_libs/sortmerna/cmph -name '*.h' -printf "%f\n"
# OpenCV only uses fraction of these files
set(CMPH_HDRS
    ${CMAKE_SOURCE_DIR}/cmph/bdz.h
	${CMAKE_SOURCE_DIR}/cmph/bdz_ph.h
	${CMAKE_SOURCE_DIR}/cmph/bdz_structs.h
	${CMAKE_SOURCE_DIR}/cmph/bdz_structs_ph.h
	${CMAKE_SOURCE_DIR}/cmph/bitbool.h
	${CMAKE_SOURCE_DIR}/cmph/bmz.h
	${CMAKE_SOURCE_DIR}/cmph/bmz8.h
	${CMAKE_SOURCE_DIR}/cmph/bmz8_structs.h
	${CMAKE_SOURCE_DIR}/cmph/bmz_structs.h
	${CMAKE_SOURCE_DIR}/cmph/brz.h
	${CMAKE_SOURCE_DIR}/cmph/brz_structs.h
	${CMAKE_SOURCE_DIR}/cmph/buffer_entry.h
	${CMAKE_SOURCE_DIR}/cmph/buffer_manager.h
	${CMAKE_SOURCE_DIR}/cmph/chd.h
	${CMAKE_SOURCE_DIR}/cmph/chd_ph.h
	${CMAKE_SOURCE_DIR}/cmph/chd_structs.h
	${CMAKE_SOURCE_DIR}/cmph/chd_structs_ph.h
	${CMAKE_SOURCE_DIR}/cmph/chm.h
	${CMAKE_SOURCE_DIR}/cmph/chm_structs.h
	${CMAKE_SOURCE_DIR}/cmph/cmph.h
	${CMAKE_SOURCE_DIR}/cmph/cmph_structs.h
	${CMAKE_SOURCE_DIR}/cmph/cmph_time.h
	${CMAKE_SOURCE_DIR}/cmph/cmph_types.h
	${CMAKE_SOURCE_DIR}/cmph/compressed_rank.h
	${CMAKE_SOURCE_DIR}/cmph/compressed_seq.h
	${CMAKE_SOURCE_DIR}/cmph/debug.h
	${CMAKE_SOURCE_DIR}/cmph/fch.h
	${CMAKE_SOURCE_DIR}/cmph/fch_buckets.h
	${CMAKE_SOURCE_DIR}/cmph/fch_structs.h
	${CMAKE_SOURCE_DIR}/cmph/graph.h
	${CMAKE_SOURCE_DIR}/cmph/hash.h
	${CMAKE_SOURCE_DIR}/cmph/hash_state.h
	${CMAKE_SOURCE_DIR}/cmph/jenkins_hash.h
	${CMAKE_SOURCE_DIR}/cmph/miller_rabin.h
	${CMAKE_SOURCE_DIR}/cmph/select.h
	${CMAKE_SOURCE_DIR}/cmph/select_lookup_tables.h
	${CMAKE_SOURCE_DIR}/cmph/vqueue.h
	${CMAKE_SOURCE_DIR}/cmph/vstack.h
)

# find C:/a03_libs/sortmerna/cmph -name '*.c' -printf "%f\n"
# OpenCV only uses fraction of these files
set(CMPH_SRCS
	${CMAKE_SOURCE_DIR}/cmph/bdz.c
	${CMAKE_SOURCE_DIR}/cmph/bdz_ph.c
	${CMAKE_SOURCE_DIR}/cmph/bmz.c
	${CMAKE_SOURCE_DIR}/cmph/bmz8.c
	${CMAKE_SOURCE_DIR}/cmph/brz.c
	${CMAKE_SOURCE_DIR}/cmph/buffer_entry.c
	${CMAKE_SOURCE_DIR}/cmph/buffer_manager.c
	${CMAKE_SOURCE_DIR}/cmph/chd.c
	${CMAKE_SOURCE_DIR}/cmph/chd_ph.c
	${CMAKE_SOURCE_DIR}/cmph/chm.c
	${CMAKE_SOURCE_DIR}/cmph/cmph.c
	${CMAKE_SOURCE_DIR}/cmph/cmph_structs.c
	${CMAKE_SOURCE_DIR}/cmph/compressed_rank.c
	${CMAKE_SOURCE_DIR}/cmph/compressed_seq.c
	${CMAKE_SOURCE_DIR}/cmph/fch.c
	${CMAKE_SOURCE_DIR}/cmph/fch_buckets.c
	${CMAKE_SOURCE_DIR}/cmph/graph.c
	${CMAKE_SOURCE_DIR}/cmph/hash.c
	${CMAKE_SOURCE_DIR}/cmph/jenkins_hash.c
	${CMAKE_SOURCE_DIR}/cmph/miller_rabin.c
	${CMAKE_SOURCE_DIR}/cmph/select.c
	${CMAKE_SOURCE_DIR}/cmph/vqueue.c
	${CMAKE_SOURCE_DIR}/cmph/vstack.c
)

#message("CMPH_SRCS= ${CMPH_SRCS}\nCMPH_PRIVATE_HDRS= ${CMPH_PRIVATE_HDRS}")

add_library(cmph OBJECT ${CMPH_SRCS} ${CMPH_HDRS})

target_include_directories(cmph PUBLIC ${CMAKE_SOURCE_DIR}/cmph)