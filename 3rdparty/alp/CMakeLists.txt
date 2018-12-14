#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 10, 2017
#
# CMake file for sortmerna/3rdparty/alp
#===============================================================================

cmake_minimum_required(VERSION 3.12)
project(alp CXX)

set(ALP_HDRS
	njn_approx.hpp
	njn_doubletype.hpp
	njn_dynprogprob.hpp
	njn_dynprogproblim.hpp
	njn_dynprogprobproto.hpp
	njn_function.hpp
	njn_integer.hpp
	njn_ioutil.hpp
	njn_localmaxstat.hpp
	njn_localmaxstatmatrix.hpp
	njn_localmaxstatutil.hpp
	njn_matrix.hpp
	njn_memutil.hpp
	njn_random.hpp
	njn_root.hpp
	njn_uniform.hpp
	njn_vector.hpp
	sls_alignment_evaluer.hpp
	sls_alp.hpp
	sls_alp_data.hpp
	sls_alp_regression.hpp
	sls_alp_sim.hpp
	sls_basic.hpp
	sls_normal_distr_array.hpp
	sls_pvalues.hpp
)

set (ALP_SRCS
	njn_dynprogprob.cpp
	njn_dynprogproblim.cpp
	njn_dynprogprobproto.cpp
	njn_ioutil.cpp
	njn_localmaxstat.cpp
	njn_localmaxstatmatrix.cpp
	njn_localmaxstatutil.cpp
	njn_random.cpp
	sls_alignment_evaluer.cpp
	sls_alp.cpp
	sls_alp_data.cpp
	sls_alp_regression.cpp
	sls_alp_sim.cpp
	sls_basic.cpp
	sls_pvalues.cpp
)


add_Library(alp OBJECT ${ALP_HDRS} ${ALP_SRCS})
target_include_directories(alp PUBLIC ${CMAKE_SOURCE_DIR}/3rdparty/alp)

#add_Library(alp INTERFACE)
#target_sources(alp INTERFACE ${ALP_HDRS})
