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


/**
 * Queue for Reads' records. Concurrently accessed by the Reader (producer) and the Processors (consumers)
 */
class ReadsQueue {
	std::string id;
	std::queue<Read> recs; // shared: Reader & Processors, Writer & Processors
	int capacity; // max size of the queue
				  //	int queueSizeAvr; // average size of the queue
	bool doneAdding; // flag indicating no more records will be added. Shared.
	int numPushed = 0; // shared
	int numPopped = 0; // shared
	int numPushers; // counter of threads that push reads on this queue. Used to calculate if the pushing is over.

	std::mutex lock;
	std::condition_variable cv;

	std::stringstream ss;

public:
	ReadsQueue(std::string id, int capacity, int numPushers)
		:
		id(id),
		capacity(capacity),
		doneAdding(false),
		numPushers(numPushers)
	{
		ss << id << " created\n";
		std::cout << ss.str(); ss.str("");
	}
	~ReadsQueue() {
		ss << "Destructor called on " << id << "  recs.size= " << recs.size() << " pushed: " << numPushed << "  popped: " << numPopped << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	void push(Read & readsrec) {
		std::unique_lock<std::mutex> queueLock(lock);
		cv.wait(queueLock, [this] {return recs.size() < capacity;});
		recs.push(std::move(readsrec));
		++numPushed;
		ss << id << " Pushed id: " << readsrec.id << " header: " << readsrec.header << " sequence: " << readsrec.sequence << std::endl;
		std::cout << ss.str(); ss.str("");
		//		l.unlock();
		cv.notify_one();
	}

	Read pop() {
		std::unique_lock<std::mutex> queueLock(lock);
		cv.wait(queueLock, [this] { return doneAdding || !recs.empty();}); //  if False - keep waiting, else - proceed.
		Read rec;
		if (!recs.empty()) {
			rec = recs.front();
			recs.pop();
			++numPopped;
			//if (numPopped % 10000 == 0)
			//{
			ss << id << " Popped id: " << rec.id << " header: " << rec.header << " sequence: " << rec.sequence << std::endl;
			std::cout << ss.str(); ss.str("");
			//				printf("Thread %s Pushed: %d Popped: %d\n", ss.str().c_str(), numPushed, numPopped);
			//printf("\rThread %s Pushed: %d Popped: %d", ss.str().c_str(), numPushed, numPopped);
			//}
		}
		//	l.unlock(); // probably redundant. The lock is released when function returns
		cv.notify_one();
		return rec;
	}

	void mDoneAdding() {
		std::lock_guard<std::mutex> l(lock);
		--numPushers;
		if (numPushers == 0)
			doneAdding = true;
		cv.notify_one(); // otherwise pop can stuck not knowing the adding stopped
	}

	bool isDone() { return doneAdding && recs.empty(); }
}; // ~class ReadsQueue
