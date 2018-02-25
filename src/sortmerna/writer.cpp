/**
* FILE: writer.cpp
* Created: Nov 26, 2017 Sun
*/
#include <iomanip>
#include <sstream>
#include <chrono>
#include <iostream>

#include "writer.hpp"


// write read alignment results to disk using e.g. RocksDB
void Writer::write()
{
	std::stringstream ss;
	ss << "Writer " << id << " thread " << std::this_thread::get_id() << " started" << std::endl;
	std::cout << ss.str(); ss.str("");

	auto t = std::chrono::high_resolution_clock::now();
	int numPopped = 0;
	for (;!writeQueue.isDone();) {
		Read read = writeQueue.pop();
		if (!read.isValid || read.isEmpty) continue;
		++numPopped;
		//std::string matchResultsStr = read.matchesToJson();
		std::string readstr = read.toString();
		if (readstr.size() > 0)
			kvdb.put(std::to_string(read.id), readstr);
	}
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
	ss << std::setprecision(2) << std::fixed << "Writer " << id << " thread " << std::this_thread::get_id() 
		<< " done. Elapsed time: " << elapsed.count() << " s Reads written: " << numPopped << std::endl;
	std::cout << ss.str(); ss.str("");
} // Writer::write