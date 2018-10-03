#pragma once

#include <string>

#include "readsqueue.hpp"
#include "kvdb.hpp"

class Writer {
public:
	Writer(std::string id, ReadsQueue & writeQueue, KeyValueDatabase & kvdb, Runopts & opts)
		: id(id), writeQueue(writeQueue), kvdb(kvdb), opts(opts) {}
	~Writer() {}

	void operator()() { write(); }
	void write();
private:
	std::string id;
	ReadsQueue & writeQueue; // shared with Processor
	KeyValueDatabase & kvdb; // key-value database path (from Options)
	Runopts & opts;
};
