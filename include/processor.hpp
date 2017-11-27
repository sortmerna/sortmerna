#pragma once
/**
* FILE: processor.hpp
* Created: Nov 06, 2017 Mon
*/
#include <string>

#include "readsqueue.hpp"
#include "readstats.hpp"
#include "output.hpp"
#include "index.hpp"
#include "references.hpp"

class Processor {
public:
	Processor(std::string id,
		ReadsQueue & readQueue,
		ReadsQueue & writeQueue,
		Readstats & readstats,
		Index & index,
		References & refs,
		Output & output,
		std::function<void(Index & index, References & refs, Output & output, Readstats & readstats, Read read)> callback
	) :
		id(id),
		readQueue(readQueue),
		writeQueue(writeQueue),
		readstats(readstats),
		index(index),
		refs(refs),
		output(output),
		callback(callback) {}

	void operator()() { process(); }
	void process(); // TODO: make private?
private:
	std::string id;
	ReadsQueue & readQueue;
	ReadsQueue & writeQueue;
	Readstats & readstats;
	References & refs;
	Output & output;
	Index & index;
	std::function<void(Index & index, References & refs, Output & output, Readstats & readstats, Read & read)> callback;
};
