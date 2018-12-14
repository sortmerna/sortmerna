#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 08, 2017
#
# CMake file for sortmerna/cmph
#===============================================================================

cmake_minimum_required(VERSION 3.12)
project(cmph CXX)

#AUX_SOURCE_DIRECTORY(${CMAKE_SOURCE_DIR}/cmph CMPH_SRCS)
#include_directories(${CMAKE_SOURCE_DIR}/cmph)              # not recommended practice. 
#file(GLOB CMPH_HDRS ${CMAKE_SOURCE_DIR}/cmph "*.h" )       # not recommended practice. List files explicitly

# find C:/a03_libs/sortmerna/cmph -name '*.h' -printf "%f\n"
# OpenCV only uses fraction of these files
set(CMPH_HDRS
    bdz.h
	bdz_ph.h
	bdz_structs.h
	bdz_structs_ph.h
	bitbool.h
	bmz.h
	bmz8.h
	bmz8_structs.h
	bmz_structs.h
	brz.h
	brz_structs.h
	buffer_entry.h
	buffer_manager.h
	chd.h
	chd_ph.h
	chd_structs.h
	chd_structs_ph.h
	chm.h
	chm_structs.h
	cmph.h
	cmph_structs.h
	cmph_time.h
	cmph_types.h
	compressed_rank.h
	compressed_seq.h
	debug.h
	fch.h
	fch_buckets.h
	fch_structs.h
	graph.h
	hash.h
	hash_state.h
	jenkins_hash.h
	miller_rabin.h
	select.h
	select_lookup_tables.h
	vqueue.h
	vstack.h
)

# find C:/a03_libs/sortmerna/cmph -name '*.c' -printf "%f\n"
# OpenCV only uses fraction of these files
set(CMPH_SRCS
	bdz.c
	bdz_ph.c
	bmz.c
	bmz8.c
	brz.c
	buffer_entry.c
	buffer_manager.c
	chd.c
	chd_ph.c
	chm.c
	cmph.c
	cmph_structs.c
	compressed_rank.c
	compressed_seq.c
	fch.c
	fch_buckets.c
	graph.c
	hash.c
	jenkins_hash.c
	miller_rabin.c
	select.c
	vqueue.c
	vstack.c
)

#message("CMPH_SRCS= ${CMPH_SRCS}\nCMPH_PRIVATE_HDRS= ${CMPH_PRIVATE_HDRS}")

add_library(cmph OBJECT ${CMPH_SRCS} ${CMPH_HDRS})
target_include_directories(cmph PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/cmph)