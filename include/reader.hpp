#pragma once
/**
 * FILE: reader.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Reads 'Reads file' and pushes the records to a shared queue for further processing
 */

#include <string>
#include <fstream> // std::ifstream

#include "readsqueue.hpp"
#include "kvdb.hpp"
#include "options.hpp"

class Read; // forward

// reads Reads and Readstats files, generates Read objects and pushes them onto ReadsQueue
class Reader {
public:
	Reader(std::string id, Runopts & opts, ReadsQueue & readQueue, KeyValueDatabase & kvdb, int loopCount)
		: 
		id(id),
		opts(opts),
		readQueue(readQueue),
		kvdb(kvdb),
		loopCount(loopCount)
	{}

	void operator()() { read(); }
	void read();
	static bool loadReadByIdx(Runopts & opts, Read & read);
	static bool loadReadById(Runopts & opts, Read & read);
private:
	std::string id;
	int loopCount; // counter of processing iterations.
	Runopts & opts;
	ReadsQueue & readQueue; // shared with Processor
	KeyValueDatabase & kvdb; // key-value database
};