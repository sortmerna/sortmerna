/*
 * file: thread_errors.hpp
 *
 * An exception that escapes a std::thread entry function calls
 * std::terminate(), which kills the whole process - including a program
 * that uses SortMeRNA as a library. ThreadErrors wraps worker entry
 * functions so that the first exception thrown by any worker is captured
 * and can be rethrown by the spawning thread after join(). The worker also
 * logs through the spawning thread's log callback (see common.hpp).
 *
 * Usage:
 *   ThreadErrors errs;
 *   for (...) tpool.emplace_back(errs.spawn(worker, arg1, std::ref(arg2)));
 *   for (auto& t : tpool) t.join();
 *   errs.rethrow(); // no-op if every worker succeeded
 */

#pragma once

#include <exception>
#include <mutex>
#include <thread>
#include <tuple>
#include <utility>

#include "common.hpp" // smr_tl_log_callback

namespace sortmerna {

class ThreadErrors {
public:
	/* Start a thread running f(args...), recording any exception it throws
	 * instead of letting it escape. Arguments are forwarded as std::thread
	 * would: use std::ref for reference parameters. */
	template<typename F, typename... Args>
	std::thread spawn(F&& f, Args&&... args)
	{
		return std::thread(
			[this, log_cb = ::smr_tl_log_callback, log_ud = ::smr_tl_log_user_data,
			 fn = std::forward<F>(f), tup = std::make_tuple(std::forward<Args>(args)...)]() mutable {
				::smr_tl_log_callback = log_cb;
				::smr_tl_log_user_data = log_ud;
				try {
					std::apply(fn, std::move(tup));
				}
				catch (...) {
					record(std::current_exception());
				}
			});
	}

	/* Rethrow the first exception recorded by a worker, if any. Call only
	 * after all workers spawned through this object have been joined. */
	void rethrow()
	{
		std::exception_ptr e;
		{
			std::lock_guard<std::mutex> lk(mtx_);
			std::swap(e, first_);
		}
		if (e) std::rethrow_exception(e);
	}

private:
	void record(std::exception_ptr e)
	{
		std::lock_guard<std::mutex> lk(mtx_);
		if (!first_) first_ = e;
	}

	std::mutex mtx_;
	std::exception_ptr first_;
};

} // namespace sortmerna
