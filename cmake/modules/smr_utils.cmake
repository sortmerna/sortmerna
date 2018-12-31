#==============================================================================
# FILE: smr_utils.cmake
# Created: Dec 31, 2018 Mon
#==============================================================================

#
# Edit rocksdb/thirdparty.inc file to configure ZLIB dependency
#
# CMake REGEX only handles 'greedy' matches making it necessary to process the file
# line by line when using 'relaxed' matches. For comparison 'rocksdb_edit_3rdparty_inc_2'
# uses stricter REGEX allowing to search the whole file in one shot.
#
# TODO: This function could be made more compact by calling a macro from 'foreach' loop.
# 
function(rocksdb_edit_3rdparty_inc)
    file(RENAME ${ROCKSDB_SRC}/thirdparty.inc ${ROCKSDB_SRC}/_thirdparty.inc)
    file(READ ${ROCKSDB_SRC}/_thirdparty.inc file_in)

    file(STRINGS ${ROCKSDB_SRC}/_thirdparty.inc lines)
    list(LENGTH lines _len)
    message("list length: ${_len}")
    list(GET lines 0 _line)
    message("line: ${_line}")

    set(zlib_home_done 0)
    set(zlib_include_done 0)
    set(zlib_lib_debug_done 0)
    set(zlib_lib_release_done 0)

    foreach(_line IN LISTS lines)
        # ERROR: string sub-command REGEX, mode MATCH needs at least 5 arguments total to command.
        # Seems to be a CMake bug. Doesn't work in the loop
        # Erroneously created issue at https://github.com/facebook/rocksdb/issues/4826 - closed
        # Opened issue: https://gitlab.kitware.com/cmake/cmake/issues/18737
        # The problem was the unquoted ${_line} i.e. always use "${_line}" to generate 
        # an empty string in case _line is Empty/Null
        if(NOT zlib_home_done)
            string(REGEX MATCH "set\\(ZLIB_HOME (.+)\\)" _mm "${_line}")
            if(_mm)
                message("mm: ${_mm}")
                message("CMAKE_MATCH_1: ${CMAKE_MATCH_1}")
                string(LENGTH "${CMAKE_MATCH_1}" _len)
                if(_len GREATER 0)
                    string(REPLACE "${CMAKE_MATCH_1}" "${ZLIB_ROOT}" _out "${_line}")
                    set(zlib_home_done 1)
                    message("ZLIB_HOME Done _out: ${_out}")
                    file(APPEND ${ROCKSDB_SRC}/thirdparty.inc "${_out}\n")
                    unset(_out)
                    unset(_mm)
                    unset(_len)
                    continue()
                endif()
            endif()
        endif()
        if(NOT zlib_include_done)
            string(REGEX MATCH "set\\(ZLIB_INCLUDE (.+)\\)" _mm "${_line}")
            if(_mm)
                string(LENGTH "${CMAKE_MATCH_1}" _len)
                if(_len GREATER 0)
                    string(REPLACE "${CMAKE_MATCH_1}" "\${ZLIB_HOME}/include" _out "${_line}")
                    set(zlib_include_done 1)
                    message("ZLIB_INCLUDE Done _out: ${_out}")
                    file(APPEND ${ROCKSDB_SRC}/thirdparty.inc "${_out}\n")
                    unset(_out)
                    unset(_mm)
                    unset(_len)
                    continue()
                endif()
            endif()
        endif()
        if(NOT zlib_lib_debug_done)
            string(REGEX MATCH "set\\(ZLIB_LIB_DEBUG (.+)\\)" _mm "${_line}")
            if(_mm)
                string(LENGTH "${CMAKE_MATCH_1}" _len)
                if(_len GREATER 0)
                    if(WIN32)
                        string(REPLACE "${CMAKE_MATCH_1}" "\${ZLIB_HOME}/lib/zlibstaticd.lib" _out "${_line}")
                    else()
                        string(REPLACE "${CMAKE_MATCH_1}" "\${ZLIB_HOME}/lib/libz.a" _out "${_line}")
                    endif()
                    set(zlib_lib_debug_done 1)
                    message("ZLIB_LIB_DEBUG Done _out: ${_out}")
                    file(APPEND ${ROCKSDB_SRC}/thirdparty.inc "${_out}\n")
                    unset(_out)
                    unset(_mm)
                    unset(_len)
                    continue()
                endif()
            endif()
        endif()
        if(NOT zlib_lib_release_done)
            string(REGEX MATCH "set\\(ZLIB_LIB_RELEASE (.+)\\)" _mm "${_line}")
            if(_mm)
                string(LENGTH "${CMAKE_MATCH_1}" _len)
                if(_len GREATER 0)
                    if(WIN32)
                        string(REPLACE "${CMAKE_MATCH_1}" "\${ZLIB_HOME}/lib/zlibstatic.lib" _out "${_line}")
                    else()
                        string(REPLACE "${CMAKE_MATCH_1}" "\${ZLIB_HOME}/lib/libz.a" _out "${_line}")
                    endif()
                    set(zlib_lib_release_done 1)
                    message("ZLIB_LIB_RELEASE Done _out: ${_out}")
                    file(APPEND ${ROCKSDB_SRC}/thirdparty.inc "${_out}\n")
                    unset(_out)
                    unset(_mm)
                    unset(_len)
                    continue()
                endif()
            endif()
        endif()
        file(APPEND ${ROCKSDB_SRC}/thirdparty.inc "${_line}\n")
    endforeach()
endfunction(rocksdb_edit_3rdparty_inc)

#
# Edit rocksdb/thirdparty.inc file to configure ZLIB dependency
#
# NOTE: Does the same as 'rocksdb_edit_3rdparty_inc' above using stricter REGEX, which simplifies
# the code but makes it less flexible in case 'thirdparty.inc' is ever modified causing
# the strict regex to fail
#
function(rocksdb_edit_3rdparty_inc_2)
    # 1. Backup the existing file
    file(RENAME ${ROCKSDB_SRC}/thirdparty.inc ${ROCKSDB_SRC}/_thirdparty.inc)
    # 2, Read the file
    file(READ ${ROCKSDB_SRC}/_thirdparty.inc file_in)
    # 3. Replace the lines
    string(REPLACE "set(ZLIB_HOME \$ENV{THIRDPARTY_HOME}/ZLIB.Library)" "set(ZLIB_HOME ${ZLIB_ROOT})" _out ${file_in})
    string(REPLACE "set(ZLIB_INCLUDE \${ZLIB_HOME}/build/native/inc/inc)" "set(ZLIB_INCLUDE ${ZLIB_ROOT}/include)" _out ${_out})
    if(WIN32)
        string(REPLACE "set(ZLIB_LIB_DEBUG \${ZLIB_HOME}/lib/native/debug/amd64/zlib.lib)" "set(ZLIB_LIB_DEBUG ${ZLIB_ROOT}/lib/zlibstaticd.lib)" _out ${_out})
        string(REPLACE "set(ZLIB_LIB_RELEASE \${ZLIB_HOME}/lib/native/retail/amd64/zlib.lib)" "set(ZLIB_LIB_RELEASE ${ZLIB_ROOT}/lib/zlibstaticd.lib)" _out ${_out})
    else()
        string(REPLACE "set(ZLIB_LIB_DEBUG \${ZLIB_HOME}/lib/native/debug/amd64/zlib.lib)" "set(ZLIB_LIB_DEBUG ${ZLIB_ROOT}/lib/libz.a)" _out ${_out})
        string(REPLACE "set(ZLIB_LIB_RELEASE \${ZLIB_HOME}/lib/native/retail/amd64/zlib.lib" "set(ZLIB_LIB_RELEASE ${ZLIB_ROOT}/lib/libz.a)" _out ${_out})
    endif()
    # Write the new file
    file(WRITE ${ROCKSDB_SRC}/thirdparty.inc ${_out})
endfunction(rocksdb_edit_3rdparty_inc_2)