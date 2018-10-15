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

	ThreadPool(int numThreads) : shutdown_(false), running_jobs(0)
	{
		std::stringstream ss;

		// Create the specified number of threads
		threads_.reserve(numThreads);
		for (int i = 0; i < numThreads; ++i)
			threads_.emplace_back(std::bind(&ThreadPool::threadEntry, this, i));

		ss << __func__ << ":" << __LINE__ << " initialized Pool with: [" << numThreads << "] threads" << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	~ThreadPool()
	{
		{
			// Unblock any threads and tell them to stop
			std::lock_guard <std::mutex> lmjq(job_queue_lock);
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
		std::lock_guard <std::mutex> lmjQ(job_queue_lock);
		jobs_.emplace(std::move(func));
		cv_jobs.notify_one();
	}

	// wait till no jobs running
	void waitAll()
	{
		std::unique_lock<std::mutex> lmjD(job_done_lock);
		//while (busy.load() != 0 && !jobs_.empty()) cv_done.wait(lmJobDone);
		cv_done.wait(lmjD, [this] { return running_jobs == 0 && jobs_.empty(); }); // works
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
				std::unique_lock <std::mutex> lockmJobQueue(job_queue_lock);

				// while no jobs and no shutdown - just keep waiting.
				while (!shutdown_.load() && jobs_.empty())
					cv_jobs.wait(lockmJobQueue); // works
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
				++running_jobs;
			}
			// mutex 'l' released here

			job(); // Do the job without holding any locks
			--running_jobs;
			cv_done.notify_one(); // whithout this the main thread hangs forever after calling 'waitAll'
			ss << __func__ << ":" << __LINE__ << " number of running_jobs= " << running_jobs << " jobs queue empty= " << jobs_.empty() << std::endl;
			std::cout << ss.str(); ss.str("");
		} // ~for
	} // ~threadEntry

	std::mutex job_queue_lock; // lock for pop/push on jobs_
	std::mutex job_done_lock; // lock for checking jobs_.empty and shutdown
	std::condition_variable cv_jobs;
	std::condition_variable cv_done;
	std::atomic_bool shutdown_;
	std::queue <std::function <void(void)>> jobs_;
	std::vector <std::thread> threads_;
	std::atomic_uint running_jobs; // counter of running jobs
}; // ~class ThreadPool

