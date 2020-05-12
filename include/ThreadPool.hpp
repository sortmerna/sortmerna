/**
* FILE: ThreadPool.hpp
* Created: 20170810 Thu
*
* https://stackoverflow.com/questions/26516683/reusing-thread-in-loop-c
*/

#pragma once

#include <iostream>
#include <sstream>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>
#include <chrono>
#include <atomic>

/**
*  all the pool threads are initially in a waiting state until jobs are available for execution.
*/
class ThreadPool
{
public:

	ThreadPool(int numThreads) : shutdown_(false), running_threads(0)
	{
		std::stringstream ss;

		// Create the specified number of threads
		threads_.reserve(numThreads);
		for (int i = 0; i < numThreads; ++i)
			threads_.emplace_back(std::bind(&ThreadPool::threadEntry, this, i));

		ss << STAMP << "initialized Pool with: [" << numThreads << "] threads" << std::endl << std::endl;
		std::cout << ss.str();
	}

	~ThreadPool()
	{
		{
			// Unblock any threads and tell them to stop
			std::lock_guard <std::mutex> lk(job_queue_mx);
			shutdown_ = true;
			cv_jobs.notify_all();
		}
		joinAll();
	}

	/** 
	 * lock the jobs queue and add a job
	 */
	void addJob(std::function <void(void)> func)
	{
		std::lock_guard <std::mutex> lk(job_queue_mx);
		jobs_.emplace(std::move(func));
		cv_jobs.notify_one();
	}

	// wait till no jobs running
	void waitAll()
	{
		std::unique_lock<std::mutex> lk(job_done_mx);
		cv_done.wait(lk, [this] { return running_threads.load() == 0 && jobs_.empty(); });
	}

	// Wait for all threads to stop
	void joinAll()
	{
		for (auto& thread : threads_)
			thread.join();
	}

protected:

	void threadEntry(int i)
	{
		std::function <void(void)> job;
		std::stringstream ss;

		for (;;)
		{
			{
				std::unique_lock <std::mutex> lk(job_queue_mx);

				// while no jobs and no shutdown - just keep waiting.
				while (!shutdown_.load() && jobs_.empty())
					cv_jobs.wait(lk); // works
				//cv_jobs.wait(jqLock, [this] { return !shutdown_.load() && jobs_.empty(); });

				if (jobs_.empty()) // only get here on shutdown = true
				{
					// No jobs to do and shutting down
					ss << "Thread  " << std::this_thread::get_id() << " job done" << std::endl;
					std::cout << ss.str(); ss.str("");
					return;
				}

				job = std::move(jobs_.front());
				jobs_.pop();
				++running_threads;
			}
			// mutex 'job_queue_mx' released here

			job(); // Do the job without holding any locks
			--running_threads;
			ss << STAMP << "number of running_threads= " << running_threads.load() << " jobs queue is empty= " << jobs_.empty() << std::endl;
			std::cout << ss.str(); ss.str("");
			cv_done.notify_one(); // wake up the main thread waiting in 'waitAll'. Keep it here for multi-reference cases.
		} // ~for
	} // ~threadEntry

	std::mutex job_queue_mx; // lock for pop/push on jobs_
	std::mutex job_done_mx; // lock for checking jobs_.empty and shutdown
	std::condition_variable cv_jobs;
	std::condition_variable cv_done;
	std::atomic_bool shutdown_;
	std::atomic_uint running_threads; // counter of running threads
	std::vector <std::thread> threads_;
	std::queue <std::function <void(void)>> jobs_;
}; // ~class ThreadPool

