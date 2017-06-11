#===============================================================================
# FILE: CMakeLists.txt
# Created: Jun 10, 2017
#
# CMake file for sortmerna/3rdparty/alp
#===============================================================================

cmake_minimum_required(VERSION 3.0)
project(alp CXX)

set(ALP_HDRS
	${CMAKE_SOURCE_DIR}/alp/njn_approx.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_doubletype.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_dynprogprob.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_dynprogprob.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_dynprogproblim.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_dynprogproblim.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_dynprogprobproto.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_dynprogprobproto.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_function.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_integer.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_ioutil.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_ioutil.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_localmaxstat.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_localmaxstat.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_localmaxstatmatrix.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_localmaxstatmatrix.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_localmaxstatutil.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_localmaxstatutil.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_matrix.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_memutil.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_random.cpp
	${CMAKE_SOURCE_DIR}/alp/njn_random.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_root.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_uniform.hpp
	${CMAKE_SOURCE_DIR}/alp/njn_vector.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_alignment_evaluer.cpp
	${CMAKE_SOURCE_DIR}/alp/sls_alignment_evaluer.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp.cpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp_data.cpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp_data.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp_regression.cpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp_regression.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp_sim.cpp
	${CMAKE_SOURCE_DIR}/alp/sls_alp_sim.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_basic.cpp
	${CMAKE_SOURCE_DIR}/alp/sls_basic.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_normal_distr_array.hpp
	${CMAKE_SOURCE_DIR}/alp/sls_pvalues.cpp
	${CMAKE_SOURCE_DIR}/alp/sls_pvalues.hpp
)

add_Library(alp INTERFACE)

target_sources(alp INTERFACE ${ALP_HDRS})

target_include_directories(alp INTERFACE 
	${CMAKE_SOURCE_DIR}/alp
)
