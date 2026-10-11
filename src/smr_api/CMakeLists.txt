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

#===============================================================================
# Self-contained archive for distribution
#
# The archive above is only usable from inside this build tree: the linker
# picks its dependencies up from the other build directories. The installed
# library has no build tree to draw on, so the dependencies are merged into a
# single archive here, leaving an application to link nothing but
#
#   cc app.c -I<prefix>/include <prefix>/lib/libsmr_api.a -lstdc++ -lpthread -ldl -lm
#
# rapidgzip is absent from the list below on purpose: rapidgzip::rapidgzip is
# an INTERFACE library whose sources are headers, so its code is compiled into
# smr_objs and is already inside libsmr_api.a.
#===============================================================================

if(NOT ROCKSDB_USE_STATIC_LIBS)
    message(FATAL_ERROR
        "smr_api needs a static RocksDB to produce a self-contained "
        "libsmr_api.a - a shared library cannot be merged into an archive. "
        "Configure with -DROCKSDB_USE_STATIC_LIBS=ON.")
endif()

foreach(lib ${SMR_ROCKSDB_ARCHIVE} ${SMR_ZLIB_ARCHIVE} ${PARASAIL_LIB})
    if(NOT lib MATCHES "\\.a$")
        message(FATAL_ERROR
            "smr_api needs static dependencies to produce a self-contained "
            "libsmr_api.a, but '${lib}' is not an archive.")
    endif()
endforeach()

set(SMR_API_BUNDLE_DIR ${CMAKE_CURRENT_BINARY_DIR}/bundle/$<CONFIG>)
set(SMR_API_BUNDLE     ${SMR_API_BUNDLE_DIR}/libsmr_api.a)
set(SMR_API_MRI        ${SMR_API_BUNDLE_DIR}/smr_api.mri)

# Archives to merge, in the order the linker sees them when it links
# test_smr_api. Order matters for one pair only: libz.a (zlib) and
# libzlibstatic.a (the ZLIB_COMPAT build of zlib-ng that rapidgzip brings in)
# define the same symbol names. An archive's symbol index resolves a name to
# the first member defining it, so zlib wins, exactly as it does today.
# Elsewhere the order is irrelevant - the linker iterates within a single
# archive until closure, which is also why the bundle needs no equivalent of
# the repeated libz.a in the link line.
set(SMR_API_BUNDLE_LIBS
    $<TARGET_FILE:smr_api>
    ${SMR_ROCKSDB_ARCHIVE}
    ${SMR_ZLIB_ARCHIVE}
)
set(SMR_API_BUNDLE_DEPS smr_api build_version alp)
if(TARGET zlibstatic)
    list(APPEND SMR_API_BUNDLE_LIBS $<TARGET_FILE:zlibstatic>)
    list(APPEND SMR_API_BUNDLE_DEPS zlibstatic)
endif()
list(APPEND SMR_API_BUNDLE_LIBS
    ${PARASAIL_LIB}
    $<TARGET_FILE:build_version>
    $<TARGET_FILE:alp>
)

# 'ar -M' reads an MRI script: create, one addlib per archive, save, end.
# addlib is used rather than extracting and re-archiving because 'ar x'
# collapses the member basenames that collide across rocksdb and parasail.
list(JOIN SMR_API_BUNDLE_LIBS "\naddlib " SMR_API_ADDLIBS)
file(GENERATE
    OUTPUT ${SMR_API_MRI}
    CONTENT "create ${SMR_API_BUNDLE}\naddlib ${SMR_API_ADDLIBS}\nsave\nend\n"
)

add_custom_command(
    OUTPUT ${SMR_API_BUNDLE}
    COMMAND ${CMAKE_COMMAND} -E make_directory ${SMR_API_BUNDLE_DIR}
    COMMAND ${CMAKE_COMMAND} -DAR=${CMAKE_AR} -DMRI=${SMR_API_MRI}
            -P ${CMAKE_SOURCE_DIR}/cmake/bundle_static_libs.cmake
    DEPENDS ${SMR_API_BUNDLE_DEPS} ${SMR_API_MRI}
            ${SMR_ROCKSDB_ARCHIVE} ${SMR_ZLIB_ARCHIVE} ${PARASAIL_LIB}
    COMMENT "Bundling libsmr_api.a with its static dependencies"
    VERBATIM
)

add_custom_target(smr_api_bundle ALL DEPENDS ${SMR_API_BUNDLE})

# tests/smr_api links the bundle to check that it is self-contained
set(SMR_API_BUNDLE "${SMR_API_BUNDLE}" PARENT_SCOPE)

# Install the bundle, not the target: $<TARGET_FILE:smr_api> is incomplete.
# Destinations are relative so that CPack stages them.
install(FILES ${SMR_API_BUNDLE}
    DESTINATION lib
    COMPONENT devel
)

# By name, not install(DIRECTORY include/): include/ is a flat mix of this one
# public header and the internal C++ headers.
install(FILES ${CMAKE_SOURCE_DIR}/include/smr_api.h
    DESTINATION include
    COMPONENT devel
)
