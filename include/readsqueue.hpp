#pragma once
/**
* FILE: readsqueue.hpp
* Created: Nov 06, 2017 Mon
*/
#include <string>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <sstream>

#include "read.hpp"
#include "concurrentqueue.h"


/**
 * Queue for Reads' records. Concurrently accessed by the Reader (producer) and the Processors (consumers)
 */
class ReadsQueue {
	std::string id;
	int capacity; // max size of the queue
	int numPushed; // shared
	int numPopped; // shared
	int numPushers; // counter of threads that push reads on this queue. Flags when the pushing is over.
	//int queueSizeAvr; // average size of the queue
#ifdef LOCKQUEUE
	std::queue<Read> recs; // shared: Reader & Processors, Writer & Processors
#else
	moodycamel::ConcurrentQueue<Read> recs; // lockless queue
#endif

	std::mutex qlock; // lock for push/pop on queue
	std::condition_variable cvQueue;

	std::stringstream ss;

public:
	ReadsQueue(std::string id, int capacity, int numPushers)
		:
		id(id),
		capacity(capacity),
		numPushed(0),
		numPopped(0),
		numPushers(numPushers),
#ifndef LOCKQUEUE
		recs(capacity) // set initial capacity
#endif
	{
		ss << id << " created" << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	~ReadsQueue() {
#ifdef LOCKQUEUE
		size_t recsize = recs.size();
#else
		size_t recsize = recs.size_approx();
#endif
		ss << "Destructor called on " << id << "  recs.size= " << recsize << " pushed: " << numPushed << "  popped: " << numPopped << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	void push(Read & rec) {
#ifdef LOCKQEUEU
		std::unique_lock<std::mutex> lmq(qlock);
		cvQueue.wait(lmq, [this] {return recs.size() < capacity;});
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

	Read pop() {
		Read rec;
#ifdef LOCKQEUEU
		std::unique_lock<std::mutex> lmq(qlock);
		cvQueue.wait(lmq, [this] { return doneAdding || !recs.empty();}); //  if False - keep waiting, else - proceed.
		if (!recs.empty()) {
			rec = recs.front();
			recs.pop();
			++numPopped;
			if (numPopped % 100000 == 0)
			{
				ss << id << " Popped id: " << rec.id << "\r"; //" Index: " << rec.lastIndex << " Part: " << rec.lastPart
			//	<< " header: " << rec.header << " sequence: " << rec.sequence << std::endl;
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

	void mDoneAdding() {
		std::lock_guard<std::mutex> lmq(qlock);
		--numPushers;
		ss << "ReadsQueue: " << id << " mDoneAdding. numPushers: " << numPushers << std::endl; std::cout << ss.str(); ss.str("");
		cvQueue.notify_one(); // otherwise pop can stuck not knowing the adding stopped
	}

	// done when no more adding and no records
	bool isDone() {
		std::lock_guard<std::mutex> lmq(qlock);
#ifdef LOCKQEUEU
		return numPushers == 0 && recs.empty();
#else
		return numPushers == 0 && recs.size_approx() == 0;
#endif
		cvQueue.notify_one(); // otherwise pop can stuck not knowing the adding stopped
	}

	// call from main thread when no other threads running
	void reset(int nPushers) {
		numPushers = nPushers;
		ss << "ReadsQueue id: " << id << " reset. numPushers: " << numPushers << std::endl; std::cout << ss.str(); ss.str("");
	}
}; // ~class ReadsQueue
