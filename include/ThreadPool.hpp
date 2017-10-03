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

	ThreadPool(int numThreads) : shutdown_(false)
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
			condVar_.notify_all();
		}

		// Wait for all threads to stop
//		std::cerr << "Joining threads" << std::endl;
		printf("Joining threads\n");
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
		condVar_.notify_one();
	}

protected:

	void threadEntry(int i)
	{
		std::function <void(void)> job;

		while (1)
		{
			{
				std::unique_lock <std::mutex> l(lock_);

				// while no jobs and no shutdown - just keep waiting.
				while (!shutdown_ && jobs_.empty())
					condVar_.wait(l);

				if (jobs_.empty())
				{
					// No jobs to do and shutting down
//					std::cerr << "Thread " << std::this_thread::get_id() << " terminates" << std::endl;
					printf("Thread %d terminates\n", std::this_thread::get_id());
					return;
				}

//				std::cerr << "Thread " << std::this_thread::get_id() << " running a job" << std::endl;
				printf("Thread %d running a job\n", std::this_thread::get_id());
				job = std::move(jobs_.front());
				jobs_.pop();
			}

			job(); // Do the job without holding any locks
		}
	} // ~threadEntry

	std::mutex lock_;
	std::condition_variable condVar_;
	bool shutdown_;
	std::queue <std::function <void(void)>> jobs_;
	std::vector <std::thread> threads_;
}; // ~class ThreadPool

