/**
* FILE: ThreadPool.hpp
* Created: 20170810 Thu
*
* https://stackoverflow.com/questions/26516683/reusing-thread-in-loop-c
*/

#pragma once

#include <iostream>
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

	ThreadPool(int numThreads) : shutdown_(false), busy(0)
	{
		// Create the specified number of threads
		threads_.reserve(numThreads);
		for (int i = 0; i < numThreads; ++i)
			threads_.emplace_back(std::bind(&ThreadPool::threadEntry, this, i));
	}

	~ThreadPool()
	{
		{
			// Unblock any threads and tell them to stop
			std::unique_lock <std::mutex> l(job_queue_lock);

			shutdown_ = true;
			cv_jobs.notify_all();
		}
		joinAll();
	}

	/** 
	 * Place a job on the queue and unblock a thread
	 */
	void addJob(std::function <void(void)> func)
	{
		std::unique_lock <std::mutex> l(job_queue_lock);
		jobs_.emplace(std::move(func));
		cv_jobs.notify_one();
	}

	// wait till no jobs running
	void ThreadPool::waitAll()
	{
		std::unique_lock<std::mutex> doneLock(job_done_lock);
		//while (busy.load() != 0 && !jobs_.empty()) cv_done.wait(doneLock);
		cv_done.wait(doneLock, [this] { return busy.load() == 0 && jobs_.empty(); }); // worked OK
		//lock.unlock();
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
				std::unique_lock <std::mutex> jqLock(job_queue_lock);

				// while no jobs and no shutdown - just keep waiting.
				while (!shutdown_.load() && jobs_.empty())	cv_jobs.wait(jqLock); // this works
				//cv_jobs.wait(jqLock, [this] { return !shutdown_ && jobs_.empty(); }); // this doesn't

				if (jobs_.empty()) // only get here on shutdown = true
				{
					// No jobs to do and shutting down
					ss << "Thread  " << std::this_thread::get_id() << " job done\n";
					std::cout << ss.str(); ss.str("");
					return;
				}

				ss << "Thread " << std::this_thread::get_id() << " running a job\n";
				std::cout << ss.str(); ss.str("");
				job = std::move(jobs_.front());
				jobs_.pop();
				++busy;
			}
			// mutex 'l' released here

			job(); // Do the job without holding any locks
			--busy;
			cv_done.notify_one(); // whithout this main thread hangs forever after calling 'waitAll'
			ss << "ThreadPool::busy= " << busy << " jobs_.empty= " << jobs_.empty() << std::endl;
			std::cout << ss.str(); ss.str("");
		} // ~for
	} // ~threadEntry

	std::mutex job_queue_lock;
	std::mutex job_done_lock;
	std::condition_variable cv_jobs;
	std::condition_variable cv_done;
	std::atomic_bool shutdown_;
	std::queue <std::function <void(void)>> jobs_;
	std::vector <std::thread> threads_;
	std::atomic_uint busy; // counter of running jobs jobs
}; // ~class ThreadPool
