/*
 * file: thread_errors.hpp
 *
 * An exception that escapes a std::thread entry function calls
 * std::terminate(), which kills the whole process - including a program
 * that uses SortMeRNA as a library. ThreadErrors wraps worker entry
 * functions so that the first exception thrown by any worker is captured
 * and can be rethrown by the spawning thread after join().
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

class ThreadErrors {
public:
	/* Start a thread running f(args...), recording any exception it throws
	 * instead of letting it escape. Arguments are forwarded as std::thread
	 * would: use std::ref for reference parameters. */
	template<typename F, typename... Args>
	std::thread spawn(F&& f, Args&&... args)
	{
		return std::thread(
			[this, fn = std::forward<F>(f), tup = std::make_tuple(std::forward<Args>(args)...)]() mutable {
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
