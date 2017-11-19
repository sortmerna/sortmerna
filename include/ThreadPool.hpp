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
		std::unique_lock<std::mutex> lock(job_done_lock);
		cv_done.wait(lock, [this] { return busy == 0; });
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
				std::unique_lock <std::mutex> l(job_queue_lock);

				// while no jobs and no shutdown - just keep waiting.
				while (!shutdown_ && jobs_.empty())
					cv_jobs.wait(l);

				if (jobs_.empty()) // only get here on shutdown = true
				{
					// No jobs to do and shutting down
//					std::cerr << "Thread " << std::this_thread::get_id() << " terminates" << std::endl;
					ss << std::this_thread::get_id();
					printf("Thread %s job done\n", ss.str().c_str());
					ss.str("");
					return;
				}

//				std::cerr << "Thread " << std::this_thread::get_id() << " running a job" << std::endl;
				ss << std::this_thread::get_id();
				printf("Thread %s running a job\n", ss.str().c_str());
				ss.str("");
				job = std::move(jobs_.front());
				jobs_.pop();
				++busy;
			}

			job(); // Do the job without holding any locks
			--busy;
			cv_done.notify_one(); // whithout this main thread calling 'waitAll' hangs forever
			//printf("ThreadPool::busy= %d jobs_.empty= %d\n", unsigned(busy), jobs_.empty());
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

