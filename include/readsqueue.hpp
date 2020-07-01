#pragma once
/**
* FILE: readsqueue.hpp
* Created: Nov 06, 2017 Mon
*/
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <sstream>
#include <atomic>

#include "common.hpp"
#include "read.hpp"

#if defined(CONCURRENTQUEUE)
#  include "concurrentqueue.h"
#elif defined(LOCKQUEUE)
#  include <queue>
#endif


/**
 * Queue for Reads' records. Concurrently accessed by the Reader (producer) and the Processors (consumers)
 */
class ReadsQueue 
{
public:
	std::string id;
	size_t capacity; // max size of the queue
	unsigned reads_tot; // total number of reads expected to be pushed/popped
	std::atomic_uint num_pushed;
	unsigned count_to_print;
	unsigned max_print_count;
	std::atomic_uint num_popped;
	std::atomic_uint num_poppers;
#if defined(CONCURRENTQUEUE)
	moodycamel::ConcurrentQueue<std::string> queue; // lockless queue
#elif defined(LOCKQUEUE)
	std::queue<Read> recs; // shared: Reader & Processors, Writer & Processors
	std::mutex qlock; // lock for push/pop on queue
	std::condition_variable cvQueue;
#endif

public:
	ReadsQueue(std::string id = "", std::size_t capacity = 100, std::size_t num_reads_tot = 0, std::size_t poppers = 0)
		:
		id(id),
		capacity(capacity),
		reads_tot(num_reads_tot),
		num_pushed(0),
		count_to_print(0),
		max_print_count(reads_tot/20),
		num_popped(0),
		num_poppers(poppers)
#ifdef CONCURRENTQUEUE
		,
		queue(capacity) // set initial capacity
#endif
	{
		INFO("created Reads queue with capacity [", capacity, "] Total reads to process: ", num_reads_tot);
	}

	//~ReadsQueue()

	/** 
	 * Synchronized. Blocks until queue has capacity for more reads
	 * pushing stops automatically upon EOF which sets is_done_push = true
	 */
	bool push(const std::string& rec) 
	{
		bool ret = false;
#if defined(CONCURRENTQUEUE)
		while (!(ret = queue.try_enqueue(rec))) {
			std::this_thread::sleep_for(std::chrono::nanoseconds(5));
		}
		if (ret) {
			num_pushed.fetch_add(1, std::memory_order_relaxed);
			++count_to_print;
		}
		if (count_to_print == max_print_count)
		{
			INFO_MEM("Thread [", std::this_thread::get_id(), "] Pushed another: ", max_print_count, ". Queue size: ", queue.size_approx());
			count_to_print = 0;

		}
		if (reads_tot == num_pushed.load(std::memory_order_relaxed))
		{
			INFO("Thread [" , std::this_thread::get_id(), "] done Push reads total: ", reads_tot, ". Queue size: ", queue.size_approx());
		}
#elif defined(LOCKQUEUE)
		std::unique_lock<std::mutex> lmq(qlock);
		cvQueue.wait(lmq, [this] { return recs.size() < capacity; });
		recs.push(std::move(rec));
		cvQueue.notify_one();
#endif
		return ret;
	}

	// synchronized
	// return false when is_done_push == true && num_pushed == num_popped
	bool pop(std::string& rec)
	{
		bool ret = false;
		//unsigned num_pop_tries = 0;
#if defined(CONCURRENTQUEUE)
		for (; !(ret = queue.try_dequeue(rec)); ) {
			std::this_thread::sleep_for(std::chrono::nanoseconds(1));
			// acquire
			// num_pushed.load(std::memory_order_relaxed) == reads_tot && queue.size_approx() == 0
			if (num_popped.load(std::memory_order_relaxed) == reads_tot) {
				break;
			}
			//INFO("Thread [", std::this_thread::get_id(), "]  Queue size: ", queue.size_approx());
			//else if (num_pop_tries > 1000) {
			//	INFO("Thread [", std::this_thread::get_id(), "] done after max Pop tries: ", num_pop_tries, ". Queue size: ", queue.size_approx());
			//	break;
			//}
		}
		if (ret)
			num_popped.fetch_add(1, std::memory_order_relaxed); // ++num_out  store  release
#elif defined(LOCKQUEUE)
		std::unique_lock<std::mutex> lmq(qlock);
		cvQueue.wait(lmq, [this] { return (pushers.load() == 0 && recs.empty()) || !recs.empty(); }); // if False - keep waiting, else - proceed.
		if (!recs.empty()) 
		{
			rec = recs.front();
			recs.pop();
			++numPopped;
			if (numPopped.load() % 100000 == 0)
			{
				std::stringstream ss;
				ss << STAMP << id << " Popped read number: " << rec.read_num << "\r";
				std::cout << ss.str();
			}
		}
		cvQueue.notify_one();
#endif
		return ret;
	}

	void reset() {
		num_pushed = 0;
		num_popped = 0;
	}

}; // ~class ReadsQueue
