#pragma once
/**
 * FILE: processor.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Callable objects designed to be run in threads
 */
#include <string>
#include <vector>
#include <functional>

// forward
class Read;
class ReadsQueue;
struct Runopts;
struct Index;
class References;
class Output;
struct Readstats;
class Refstats;

/* 
 * performs alignment
 */
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
}; // ~class Processor

/* performs post-alignment tasks like calculating statistics */
class PostProcessor {
public:
	PostProcessor(
		std::string id,
		ReadsQueue & readQueue,
		ReadsQueue & writeQueue,
		Runopts & opts,
		References & refs,
		Readstats & readstats,
		void(*callback)(Read & read, Readstats & readstats, References & refs, Runopts & opts)
	) :
		id(id),
		readQueue(readQueue),
		writeQueue(writeQueue),
		opts(opts),
		refs(refs),
		readstats(readstats),
		callback(callback)
	{}

	void operator()() { run(); }

protected:
	void run();
	void(*callback)(Read & read, Readstats & readstats, References & refs, Runopts & opts);

protected:
	std::string id;
	// callback parameters. TODO: a better way of binding. (std::bind doesn't look better)
	ReadsQueue & readQueue;
	ReadsQueue & writeQueue;
	Runopts & opts;
	References & refs;
	Readstats & readstats;
}; // ~class PostProcessor

/* generates output after alignment and post-processing are done */
class ReportProcessor {
public:
	ReportProcessor(
		std::string id,
		ReadsQueue & readQueue,
		Runopts & opts,
		References & refs, 
		Output & output, 
		Refstats & refstats,
		void(*callback)(std::vector<Read> & reads, Runopts & opts, References & refs, Refstats & refstats, Output & output)
	) :
		id(id),
		readQueue(readQueue),
		opts(opts),
		refs(refs),
		output(output),
		refstats(refstats),
		callback(callback)
	{}

	void operator()() { run(); }
	//using Processor::operator(); // otherwise explicitely generated () operator overrides superclass () operator
	//using Processor::process;
protected:
	void run();
	void(*callback)(std::vector<Read> & reads, Runopts & opts, References & refs, Refstats & refstats, Output & output);

protected:
	std::string id;
	ReadsQueue & readQueue;
	// callback parameters. TODO: a better way of binding. (std::bind doesn't look better)
	Runopts & opts; 
	References & refs;
	Refstats & refstats;
	Output & output;
}; // ~class ReportProcessor
