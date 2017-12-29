#pragma once
/**
* FILE: processor.hpp
* Created: Nov 06, 2017 Mon
*/
#include <string>
#include <functional>

//#include "processor.hpp"
//#include "readstats.hpp"
//#include "refstats.hpp"
//#include "output.hpp"

// forward
class Read;
class ReadsQueue;
struct Runopts;
struct Index;
class References;
class Output;
struct Readstats;
class Refstats;

class Processor {
public:
	Processor(
		std::string id,
		ReadsQueue & readQueue,
		ReadsQueue & writeQueue,
		Runopts & opts, 
		Index & index, 
		References & refs, 
		Output & output, 
		Readstats & readstats, 
		Refstats & refstats,
		//std::function<void(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read)> callback
		void(*callback)(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read)
	) :
		id(id),
		readQueue(readQueue),
		writeQueue(writeQueue),
		opts(opts),
		index(index),
		refs(refs),
		output(output),
		readstats(readstats),
		refstats(refstats),
		callback(callback) 
	{}

	void operator()() { run(); }

protected:
	void run();
	//std::function<void(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read)> callback;
	void(*callback)(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read);

protected:
	std::string id;
	ReadsQueue & readQueue;
	ReadsQueue & writeQueue;
	Runopts & opts; 
	Index & index; 
	References & refs; 
	Output & output; 
	Readstats & readstats; 
	Refstats & refstats;
};

class ReportProcessor {
public:
	ReportProcessor(
		std::string id,
		ReadsQueue & readQueue,
		Runopts & opts, 
		Index & index, 
		References & refs, 
		Output & output, 
		Readstats & readstats, 
		Refstats & refstats,
		void(*callback)(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read)
	) :
		id(id),
		readQueue(readQueue),
		opts(opts),
		index(index),
		refs(refs),
		output(output),
		readstats(readstats),
		refstats(refstats),
		callback(callback)
	{}

	void operator()() { run(); }
	//using Processor::operator(); // otherwise explicitely generated () operator overrides superclass () operator
	//using Processor::process;
protected:
	void run();
	void(*callback)(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read);

protected:
	std::string id;
	// callback parameters. TODO: a better way of binding. (std::bind doesn't look better)
	ReadsQueue & readQueue;
	Runopts & opts; 
	Index & index;
	References & refs;
	Output & output;
	Readstats & readstats;
	Refstats & refstats;
};
