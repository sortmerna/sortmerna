/**
* FILE: writer.cpp
* Created: Nov 26, 2017 Sun
*/
#include "writer.hpp"


// write read alignment results to disk using e.g. RocksDB
void Writer::write()
{
	//std::stringstream ss;
	//ss << std::this_thread::get_id();
	//printf("Writer thread %s started\n", ss.str().c_str());
	//ss.str("");
	std::cout << "Writer thread " << std::this_thread::get_id() << " started\n";
	auto t = std::chrono::high_resolution_clock::now();
	int numPopped = 0;
	for (;;) {
		Read read = writeQueue.pop();
		if (!read.isValid || read.isEmpty)
		{
			if (writeQueue.isDone()) break;
			else continue;
		}
		++numPopped;
		//std::string matchResultsStr = read.matchesToJson();
		std::string readstr = read.toString();
		if (readstr.size() > 0)
			kvdb.put(std::to_string(read.id), readstr);
		//if (writeQueue.isDone()) break;
	}
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
	//ss << std::this_thread::get_id();
	std::cout << std::setprecision(2) << std::fixed;
	std::cout << "Writer thread " << std::this_thread::get_id() << " done. Elapsed time: " << elapsed.count() << " s Reads written: " << numPopped << std::endl;
	//printf("Writer thread %s done. Elapsed time: %.2f s Reads written: %d\n", ss.str().c_str(), elapsed.count(), numPopped);
} // Writer::write