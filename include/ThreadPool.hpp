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

	ThreadPool(int numThreads) : shutdown_(false), busy()
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
			std::unique_lock <std::mutex> l(lock_);

			shutdown_ = true;
			cv_jobs.notify_all();
		}

		// Wait for all threads to stop
//		std::cerr << "Joining threads" << std::endl;
//		printf("Joining threads\n");
		for (auto& thread : threads_)
			thread.join();
	}

	/** 
	 * Place a job on the queue and unblock a thread
	 */
	void addJob(std::function <void(void)> func)
	{
		std::unique_lock <std::mutex> l(lock_);
		jobs_.emplace(std::move(func));
		cv_jobs.notify_one();
	}

	void ThreadPool::waitDone()
	{
		std::unique_lock<std::mutex> lock(lock_);
		// wait till no more jobs and no workers running
		cv_done.wait(lock, [this]() { return jobs_.empty() && (busy == 0); });
	}

protected:

	void threadEntry(int i)
	{
		std::function <void(void)> job;
		std::stringstream ss;

		while (1)
		{
			{
				std::unique_lock <std::mutex> l(lock_);

				// while no jobs and no shutdown - just keep waiting.
				while (!shutdown_ && jobs_.empty())
					cv_jobs.wait(l);

				if (jobs_.empty())
				{
					// No jobs to do and shutting down
//					std::cerr << "Thread " << std::this_thread::get_id() << " terminates" << std::endl;
					ss << std::this_thread::get_id();
					printf("Thread %s terminates\n", ss.str().c_str());
					ss.str("");
					--busy;
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
		}
	} // ~threadEntry

	std::mutex lock_;
	std::condition_variable cv_jobs;
	std::condition_variable cv_done;
	bool shutdown_;
	std::queue <std::function <void(void)>> jobs_;
	std::vector <std::thread> threads_;
	unsigned int busy;
}; // ~class ThreadPool

