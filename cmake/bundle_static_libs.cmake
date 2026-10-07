#===============================================================================
# file: bundle_static_libs.cmake
#
# Merge several static archives into one, by feeding an MRI script to 'ar -M'.
#
# This runs as 'cmake -P' rather than directly from add_custom_command because
# 'ar -M' reads its script from stdin and custom commands cannot redirect it.
# execute_process() can, via INPUT_FILE.
#
# Invoked as:
#   cmake -DAR=<ar> -DMRI=<script> -P bundle_static_libs.cmake
#===============================================================================

if(NOT AR)
	message(FATAL_ERROR "AR not set")
endif()
if(NOT MRI)
	message(FATAL_ERROR "MRI not set")
endif()
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
