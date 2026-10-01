#===============================================================================
# CMakeLists.txt for the smr_api static library (C library API)
#===============================================================================

add_library(smr_api STATIC smr_api.cpp)

target_include_directories(smr_api
    PUBLIC
        ${CMAKE_SOURCE_DIR}/include
)

# libsmr_api.a archives smr_api.o together with the smr_objs objects. The
# libraries below are not in the archive; linking the smr_api target pulls
# them in.
target_link_libraries(smr_api
    PUBLIC
        smr_objs
        build_version
        alp
        ${CMAKE_DL_LIBS}
)
