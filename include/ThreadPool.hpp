/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

/**
* file: ThreadPool.hpp
* created: 20170810 Thu
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

#include "common.hpp"

/**
*  all the pool threads are initially in a waiting state until jobs are available for execution.
*/
class ThreadPool
{
protected:
	std::mutex job_queue_mx; // lock for pop/push on jobs_
	std::mutex job_done_mx; // lock for checking jobs_.empty and shutdown
	std::condition_variable cv_jobs;
	std::condition_variable cv_done;
	std::atomic_bool shutdown_;
	std::atomic_uint running_threads; // counter of running threads
	std::vector <std::thread> threads_;
	std::queue <std::function <void(void)>> jobs_;

public:
	ThreadPool(int numThreads) : shutdown_(false), running_threads(0)
	{
		// Create the specified number of threads
		threads_.reserve(numThreads);
		for (int i = 0; i < numThreads; ++i) {
			threads_.emplace_back(std::bind(&ThreadPool::threadEntry, this, i));
		}

		INFO("initialized Pool with [", numThreads, "] threads\n");
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
		PRN_MEM("ThreadPool destructor done.");
	}

protected:
		void threadEntry(int i)
		{
			std::function<void(void)> job;

			for (;;)
			{
				{
					std::unique_lock<std::mutex> lk(job_queue_mx);

					// Sleep until there is a job to execute or a shutdown flag set
					while (!shutdown_.load() && jobs_.empty())
						cv_jobs.wait(lk); // works
					//cv_jobs.wait(lk, [this] { return !shutdown_.load() && jobs_.empty(); });

					if (jobs_.empty()) // No jobs to do and shutting down
					{
						INFO("Thread  ", std::this_thread::get_id(), " job done");
						return;
					}

					job = std::move(jobs_.front());
					jobs_.pop();
					++running_threads;
				}
				// mutex 'job_queue_mx' released here

				job(); // Do the job without holding any locks
				--running_threads;

				INFO("number of running_threads= ", running_threads.load(), " jobs queue is empty= ", jobs_.empty());

				cv_done.notify_one(); // wake up the main thread waiting in 'waitAll'. Keep it here for multi-reference cases.
			} // ~for
		} // ~threadEntry

public:
	/** 
	 * lock the jobs queue and add a job
	 */
	void addJob(std::function <void()> func)
	{
		std::lock_guard<std::mutex> lk(job_queue_mx);
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
}; // ~class ThreadPool

