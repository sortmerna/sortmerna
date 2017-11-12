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

public:
	ReadsQueue(std::string id, int capacity, int numPushers)
		:
		id(id),
		capacity(capacity),
		doneAdding(false),
		numPushers(numPushers)
	{
		printf("%s created\n", id.c_str());
	}
	~ReadsQueue() {
		printf("Destructor called on %s  recs.size= %zu pushed: %d  popped: %d\n",
			id.c_str(), recs.size(), numPushed, numPopped);
	}

	void push(Read & readsrec) {
		std::unique_lock<std::mutex> l(lock);
		cv.wait(l, [this] {return recs.size() < capacity;});
		recs.push(std::move(readsrec));
		++numPushed;
		printf("%s Pushed id: %d header: %s sequence: %s\n", id.c_str(), readsrec.id, readsrec.header.c_str(), readsrec.sequence.c_str());
		//		l.unlock();
		cv.notify_one();
	}

	Read pop() {
		std::unique_lock<std::mutex> l(lock);
		cv.wait(l, [this] { return doneAdding || !recs.empty();}); //  if False - keep waiting, else - proceed.
		Read rec;
		if (!recs.empty()) {
			// printf("%d Recs.size: %d\n", id, recs.size());
			rec = recs.front();
			recs.pop();
			++numPopped;
			//if (numPopped % 10000 == 0)
			//{
			printf("%s Popped id: %d header: %s sequence: %s\n", id.c_str(), rec.id, rec.header.c_str(), rec.sequence.c_str());
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
