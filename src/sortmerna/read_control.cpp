/**
 * FILE: read_control.cpp 
 * Created: Mar 12, 2019 Tue
 *
 * @copyright 2016-19 Clarity Genomics BVBA
 */
#include <string>
#include <iostream> // std::cout,cerr
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <ios> // std::ios_base
#include <chrono> // std::chrono
#include <thread>
#include <iomanip> // std::precision

#include "common.hpp"
#include "read_control.hpp"
#include "read.hpp"


ReadControl::ReadControl(Runopts & opts, ReadsQueue & readQueue, KeyValueDatabase & kvdb)
	:
	opts(opts),
	readQueue(readQueue),
	kvdb(kvdb),
	vreader()
{}

ReadControl::~ReadControl(){}

void ReadControl::run()
{
	std::stringstream ss;

	// open reads files for reading
	std::ifstream ifs(opts.readsfile, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) {
		ERR("failed to open " + opts.readsfile);
		exit(EXIT_FAILURE);
	}

	std::string rid("reader_" + std::to_string(1));
	Reader reader_fwd(rid, opts.is_gz);
	vreader.push_back(reader_fwd);
	if (opts.paired)
	{
		rid = "reader_" + std::to_string(2);
		Reader reader_rev(rid, opts.is_gz);
		vreader.push_back(reader_rev);
	}

	ss << STAMP << " thread: " << std::this_thread::get_id() << " started" << std::endl;
	std::cout << ss.str(); ss.str("");
	auto t = std::chrono::high_resolution_clock::now();

	Read read;
	bool done_fwd = false;
	bool done_rev = false;
	// loop calling Readers
	for (; !reader_fwd.is_done || (opts.paired && !vreader[1].is_done);)
	{
		if (!reader_fwd.is_done)
			read = reader_fwd.nextread(opts, ifs, kvdb);
		if (opts.paired && !vreader[1].is_done)
			read = vreader[1].nextread(opts, ifs, kvdb);

		if (!read.isEmpty)
			readQueue.push(read);
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
	readQueue.decrPushers(); // signal the reader done adding
	readQueue.notify(); // notify processor that might be waiting to pop

	ss << STAMP << " thread: " << std::this_thread::get_id() << " done. Elapsed time: "
		<< std::setprecision(2) << std::fixed << elapsed.count() << " sec Reads added: " << read.id + 1
		<< " readQueue.size: " << readQueue.size() << std::endl;
	std::cout << ss.str(); ss.str("");

} // ~ReadControl::run

// ~read_control.cpp
