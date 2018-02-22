#pragma once
/**
 * FILE: reader.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Reads 'Reads file' and pushes the records to a shared queue for further processing
 */

#include <string>
#include <fstream> // std::ifstream

#include "zlib.h"

#include "readsqueue.hpp"
#include "kvdb.hpp"
#include "options.hpp"

#define OUT_CHUNK 32768U /* out buffer size */
#define CHUNK 16384      /* file input buffer size */
#define RL_OK 0
#define RL_END 1
#define RL_ERR -1

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
	void initZstream();
	int getline(std::ifstream & ifs, std::string & line);
private:
	std::string id;
	int loopCount; // counter of processing iterations.
	Runopts & opts;
	ReadsQueue & readQueue; // shared with Processor
	KeyValueDatabase & kvdb; // key-value database

	// zlib related
	char* line_start;
	z_stream * pstrm;
	std::vector<unsigned char> z_in;
	std::vector<unsigned char> z_out;

private:
	int inf(std::ifstream & ifs);
};