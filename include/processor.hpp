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
		//std::function<void(Index & index, References & refs, Output & output, Readstats & readstats, Read & read)> callback
		void(*callback)(Index & index, References & refs, Output & output, Readstats & readstats, Read & read)
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
	virtual void process(); // TODO: make private?

protected:
	std::string id;
	ReadsQueue & readQueue;
	ReadsQueue & writeQueue;
	Readstats & readstats;
	References & refs;
	Output & output;
	Index & index;
	//std::function<void(Index & index, References & refs, Output & output, Readstats & readstats, Read & read)> callback;
	void (*callback)(Index & index, References & refs, Output & output, Readstats & readstats, Read & read);
};

class ReportProcessor:Processor {
public:
	ReportProcessor(
		std::string id,
		ReadsQueue & readQueue,
		ReadsQueue & writeQueue,
		Readstats & readstats,
		Index & index,
		References & refs,
		Output & output,
		void(*callback)(Index & index, References & refs, Output & output, Readstats & readstats, Read & read)
		//std::function<void(Index & index, References & refs, Output & output, Readstats & readstats, Read & read)> callback
	)
		: Processor(id, readQueue, writeQueue, readstats, index, refs, output, callback) {}

	using Processor::operator(); // otherwise explicitely generated () operator overrides superclass () operator
	void process();
};
