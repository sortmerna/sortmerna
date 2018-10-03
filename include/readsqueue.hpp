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

#ifdef LOCKQEUEU
#include <queue>
#else
#include "concurrentqueue.h"
#endif


/**
 * Queue for Reads' records. Concurrently accessed by the Reader (producer) and the Processors (consumers)
 */
class ReadsQueue 
{
	std::string id;
	int capacity; // max size of the queue
	std::atomic<int> numPushed; // shared
	std::atomic<int> numPopped; // shared
#ifdef LOCKQEUEU
	std::queue<Read> recs; // shared: Reader & Processors, Writer & Processors
#else
	moodycamel::ConcurrentQueue<Read> recs; // lockless queue
#endif

	std::mutex qlock; // lock for push/pop on queue
	std::condition_variable cvQueue;

	std::stringstream ss;

public:
	std::atomic_uint numPushers; // counter of threads that push reads on this queue. When zero - the pushing is over.

	ReadsQueue(std::string id, int capacity, int numPushers)
		:
		id(id),
		capacity(capacity),
		numPushed(0),
		numPopped(0),
		numPushers(numPushers)
#ifndef LOCKQEUEU
		,
		recs(capacity) // set initial capacity
#endif
	{
		ss << id << " created" << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	~ReadsQueue() {
#ifdef LOCKQEUEU
		size_t recsize = recs.size();
#else
		size_t recsize = recs.size_approx();
#endif
		ss << "Destructor called on " << id << "  recs.size= " << recsize << " pushed: " << numPushed << "  popped: " << numPopped << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	void push(Read & rec) 
	{
#ifdef LOCKQEUEU
		std::unique_lock<std::mutex> lmq(qlock);
		cvQueue.wait(lmq, [this] { return recs.size() < capacity; });
		recs.push(std::move(rec));

		//ss << id << " Pushed id: " << rec.id << " Index: " << rec.lastIndex << " Part: " << rec.lastPart 
		//	<< " header: " << rec.header << " sequence: " << rec.sequence << std::endl;
		//std::cout << ss.str(); ss.str("");

		cvQueue.notify_one();
#else
		recs.enqueue(rec);
#endif
		++numPushed;
	}

	Read pop() 
	{
		Read rec;
#ifdef LOCKQEUEU
		std::unique_lock<std::mutex> lmq(qlock);
		cvQueue.wait(lmq, [this] { return (numPushers == 0 && recs.empty()) || !recs.empty(); }); // if False - keep waiting, else - proceed.
		if (!recs.empty()) 
		{
			rec = recs.front();
			recs.pop();
			++numPopped;
			if (numPopped % 100000 == 0)
			{
				ss << id << " Popped id: " << rec.id << "\r";
				std::cout << ss.str(); ss.str("");
			}
		}
		cvQueue.notify_one();
#else
		bool found = recs.try_dequeue(rec);
		if (found) ++numPopped;
#endif
		return rec;
	}

	// done when no more adding and no records
	// TODO: not used
	bool isDone() {
#ifdef LOCKQEUEU
		std::lock_guard<std::mutex> lmq(qlock);
		bool done = (numPushers == 0 && recs.empty());
		cvQueue.notify_one(); // otherwise pop can stuck not knowing the adding stopped
#else
		bool done = (numPushers == 0 && recs.size_approx() == 0);
#endif
		return done;
	}

	// call from main thread when no other threads running
	void reset(int nPushers) {
		numPushers = nPushers;
		ss << "ReadsQueue id: " << id << " reset. numPushers: " << numPushers << std::endl; std::cout << ss.str(); ss.str("");
	}

	// TODO: not used
	size_t size()
	{
#ifdef LOCKQEUEU
		return recs.size();
#else
		return recs.size_approx();
#endif
	}

	/* 
	 * used by the Pushers to notify in cases when no Reads were ever pushed to the Write Queue 
	 */
	void notify()
	{
		cvQueue.notify_all();
	}
}; // ~class ReadsQueue
