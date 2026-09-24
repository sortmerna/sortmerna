/*
 * file: smr_log.cpp
 *
 * Definitions of the thread-local log callback declared in common.hpp and
 * used by the INFO/WARN/ERR macros. The C API (smr_api) sets them on the
 * calling thread for the duration of a call; in the sortmerna executable
 * they stay null and the macros write to stdout/stderr.
 */
#include "common.hpp"

thread_local smr_log_fn smr_tl_log_callback = nullptr;
thread_local void* smr_tl_log_user_data = nullptr;
