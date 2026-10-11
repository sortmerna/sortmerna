#===============================================================================
# file: bundle_static_libs.cmake
#
# Merge several static archives into one. The merge tool depends on the
# platform:
#
#   GNU/Linux: feed an MRI script to 'ar -M'
#     cmake -DAR=<ar> -DMRI=<script> -P bundle_static_libs.cmake
#
#   macOS: BSD ar has no MRI script mode, so merge with 'libtool -static'
#     cmake -DLIBTOOL=<libtool> -DOUT=<archive> -DLIBS=<list file> \
#           -P bundle_static_libs.cmake
#
# LIBS names the input archives one per line. Both input files are
# materialized with file(GENERATE) because the archive paths carry generator
# expressions.
#
# This runs as 'cmake -P' rather than directly from add_custom_command because
# 'ar -M' reads its script from stdin and custom commands cannot redirect it.
# execute_process() can, via INPUT_FILE.
#===============================================================================

if(AR AND MRI)
	if(NOT EXISTS "${MRI}")
		message(FATAL_ERROR "MRI script not found: ${MRI}")
	endif()

	execute_process(
		COMMAND ${AR} -M
		INPUT_FILE ${MRI}
		RESULT_VARIABLE rcode
	)

	if(NOT rcode EQUAL 0)
		message(FATAL_ERROR "'${AR} -M' failed with code ${rcode} on script ${MRI}")
	endif()
elseif(LIBTOOL AND OUT AND LIBS)
	if(NOT EXISTS "${LIBS}")
		message(FATAL_ERROR "library list not found: ${LIBS}")
	endif()

	file(STRINGS ${LIBS} libs)

	execute_process(
		COMMAND ${LIBTOOL} -static -o ${OUT} ${libs}
		RESULT_VARIABLE rcode
	)

	if(NOT rcode EQUAL 0)
		message(FATAL_ERROR "'${LIBTOOL} -static' failed with code ${rcode} on archives in ${LIBS}")
	endif()
else()
	message(FATAL_ERROR "set either AR and MRI, or LIBTOOL, OUT and LIBS")
endif()
