#===============================================================================
# CMakeLists.txt for smr_api tests
#===============================================================================

# Compile as C to verify smr_api.h is C-linkable
add_executable(test_smr_api test_smoke.c)

target_include_directories(test_smr_api
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}
)

find_package(Threads REQUIRED)

target_link_libraries(test_smr_api
    smr_api
    Threads::Threads
)

target_compile_definitions(test_smr_api
    PRIVATE
        SMR_GOLDEN_DIR="${CMAKE_SOURCE_DIR}/tests/golden"
        SMR_DATA_DIR="${CMAKE_SOURCE_DIR}/test_data"
)

add_test(NAME smr_api_smoke COMMAND test_smr_api)

#===============================================================================
# Link check for the archive that gets installed.
#
# This deliberately links the bundled libsmr_api.a as a plain file plus system
# libraries only - linking the smr_api target instead would re-add the
# dependencies from the build tree and defeat the check.
#===============================================================================

add_executable(test_smr_api_link test_link.c)

target_include_directories(test_smr_api_link
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(test_smr_api_link
    ${SMR_API_BUNDLE}
    Threads::Threads
    ${CMAKE_DL_LIBS}
    m
)

# the archive is C++; a C-only target would not pull the C++ runtime in
set_target_properties(test_smr_api_link PROPERTIES LINKER_LANGUAGE CXX)

add_dependencies(test_smr_api_link smr_api_bundle)

add_test(NAME smr_api_bundle_link COMMAND test_smr_api_link)
