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
	std::atomic_bool is_done_push; // indicates pushers done adding items to the queue
	unsigned num_reads_tot; // total number of reads expected to be put/consumed
	std::atomic_uint num_pushed; // shared
	std::atomic_uint num_popped; // shared
#if defined(CONCURRENTQUEUE)
	moodycamel::ConcurrentQueue<std::string> queue; // lockless queue
#elif defined(LOCKQUEUE)
	std::queue<Read> recs; // shared: Reader & Processors, Writer & Processors
	std::mutex qlock; // lock for push/pop on queue
	std::condition_variable cvQueue;
#endif

public:
	ReadsQueue(std::string id = "", std::size_t capacity = 100, std::size_t num_reads_tot = 0)
		:
		id(id),
		capacity(capacity),
		is_done_push(false),
		num_pushed(0),
		num_popped(0),
		num_reads_tot(num_reads_tot)
#ifdef CONCURRENTQUEUE
		,
		queue(capacity) // set initial capacity
#endif
	{
		std::stringstream ss;
		ss << STAMP << "created Reads queue with capacity [" << capacity << "]" << std::endl;
		std::cout << ss.str();
	}

	~ReadsQueue() {
		std::stringstream ss;
		ss << STAMP << "Destructor called on Reads queue. Reads added: " << num_pushed << " Reads consumed: " << num_popped << std::endl;
		std::cout << ss.str();
	}

	/** 
	 * Synchronized. Blocks until queue has capacity for more reads
	 * pushing stops automatically upon EOF which sets is_done_push = true
	 */
	bool push(std::string& rec) 
	{
#if defined(CONCURRENTQUEUE)
		while (!queue.try_enqueue(rec)) {
			std::this_thread::sleep_for(std::chrono::nanoseconds(5));
		}
		num_pushed.fetch_add(1, std::memory_order_release);
#elif defined(LOCKQUEUE)
		std::unique_lock<std::mutex> lmq(qlock);
		cvQueue.wait(lmq, [this] { return recs.size() < capacity; });
		recs.push(std::move(rec));
		cvQueue.notify_one();
#endif
		return true;
	}

	// synchronized
	// return false when is_done_push == true && num_pushed == num_popped
	bool pop(std::string& rec)
	{
		bool ret = false;
#if defined(CONCURRENTQUEUE)
		while (!(ret = queue.try_dequeue(rec))) {
			std::this_thread::sleep_for(std::chrono::nanoseconds(1));
			if (is_done_push.load(std::memory_order_acquire) && num_pushed.load(std::memory_order_acquire) == num_popped.load(std::memory_order_acquire)) {
				break;
			}
		}
		if (ret)
			num_popped.fetch_add(1, std::memory_order_release); // ++num_out  store
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
		is_done_push = false;
		num_pushed = 0;
		num_popped = 0;
	}

}; // ~class ReadsQueue
