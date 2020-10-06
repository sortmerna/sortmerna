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

//#include "concurrentqueue.h"

#include "kvdb.hpp"

// forward
class Read;
class Readfeed;
struct Runopts;
struct Index;
class References;
class Output;
struct Readstats;
class Refstats;

void align(Runopts& opts, Readstats& readstats, Output& output, Index& index, KeyValueDatabase& kvdb);
void align2(int id, Readfeed& readfeed, Index& index, References& refs, Readstats& readstats, Refstats& refstats, KeyValueDatabase& kvdb, Runopts& opts);
void postProcess(Runopts& opts, Readstats& readstats, Output& output, KeyValueDatabase& kvdb);
void postProcess2(int id, Readfeed& readfeed, Runopts& opts, References& refs, Readstats& readstats, Refstats& refstats, KeyValueDatabase& kvdb);
void computeStats(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts);
void generateReports(Runopts& opts, Readstats& readstats, Output& output, KeyValueDatabase& kvdb);
void reportsJob(Readfeed& readfeed, Runopts& opts, References& refs, Refstats& refstats, Output& output, KeyValueDatabase& kvdb);

/* 
 * performs alignment
 */
class Processor {
public:
	Processor(
		int id,
		Readfeed& readfeed,
		Runopts& opts, 
		Index& index, 
		References& refs, 
		Readstats& readstats, 
		Refstats& refstats,
		KeyValueDatabase& kvdb,
		//std::function<void(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read)> callback
		void(*callback)(Runopts& opts, Index& index, References& refs, Readstats& readstats, Refstats& refstats, Read& read, bool isLastStrand)
	) :
		id(id),
		readfeed(readfeed),
		opts(opts),
		index(index),
		refs(refs),
		readstats(readstats),
		refstats(refstats),
		kvdb(kvdb),
		callback(callback) 
	{}

	//void operator()() { run(); }

protected:
	void run();
	//std::function<void(Runopts & opts, Index & index, References & refs, Output & output, Readstats & readstats, Refstats & refstats, Read & read)> callback;
	void(*callback)(Runopts& opts, Index& index, References& refs, Readstats& readstats, Refstats& refstats, Read& read, bool isLastStrand);

protected:
	int id;
	Readfeed& readfeed;
	Runopts& opts; 
	Index& index; 
	References& refs; 
	Readstats& readstats; 
	Refstats& refstats;
	KeyValueDatabase& kvdb;
}; // ~class Processor

/* performs post-alignment tasks like calculating statistics */
class PostProcessor {
public:
	PostProcessor(
		int id,
		Readfeed& readfeed,
		Runopts& opts,
		References& refs,
		Readstats& readstats,
		Refstats& refstats,
		KeyValueDatabase& kvdb,
		void(*callback)(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts)
	) :
		id(id),
		readfeed(readfeed),
		opts(opts),
		refs(refs),
		readstats(readstats),
		refstats(refstats),
		kvdb(kvdb),
		callback(callback)
	{}

	//void operator()() { run(); }

protected:
	void run();
	void(*callback)(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts);

protected:
	int id;
	// callback parameters. TODO: a better way of binding. (std::bind doesn't look better)
	Readfeed& readfeed;
	Runopts& opts;
	References& refs;
	Readstats& readstats;
	Refstats& refstats;
	KeyValueDatabase& kvdb;
}; // ~class PostProcessor

/* generates output after alignment and post-processing are done */
class ReportProcessor {
public:
	ReportProcessor(
		int id,
		Readfeed& readfeed,
		Runopts& opts,
		References& refs, 
		Output& output, 
		Refstats& refstats,
		KeyValueDatabase& kvdb
		//void(*callback)(std::vector<Read>& reads, Runopts& opts, References& refs, Refstats& refstats, Output& output)
	) :
		id(id),
		readfeed(readfeed),
		opts(opts),
		refs(refs),
		output(output),
		refstats(refstats),
		kvdb(kvdb)
		//callback(callback)
	{}

	//void operator()() { run(); }

protected:
	void run();
	//void(*callback)(std::vector<Read> & reads, Runopts & opts, References & refs, Refstats & refstats, Output & output);
	void job(std::vector<Read>& reads, Runopts& opts, References& refs, Refstats& refstats, Output& output); // report job

protected:
	int id;
	Readfeed& readfeed;
	// callback parameters. TODO: a better way of binding. (std::bind doesn't look better)
	Runopts& opts; 
	References& refs;
	Refstats& refstats;
	Output& output;
	KeyValueDatabase& kvdb;
}; // ~class ReportProcessor
