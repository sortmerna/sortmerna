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
	bool is_paired = (opts.readfiles.size() == 2); // reads are paired i.e. 2 read files are supplied

	// open reads files for reading
	auto fwd_file = opts.readfiles[0];
	std::ifstream ifs_fwd(fwd_file, std::ios_base::in | std::ios_base::binary);

	if (!ifs_fwd.is_open()) 
	{
		ERR("failed to open file: [" + fwd_file + "]");
		exit(EXIT_FAILURE);
	}

	std::ifstream ifs_rev;
	if (opts.readfiles.size() == 2)
	{
		auto rev_file = opts.readfiles[1];
		ifs_rev.open(rev_file, std::ios_base::in | std::ios_base::binary);

		if (!ifs_rev.is_open()) {
			ERR("failed to open file: [" + rev_file + "]");
			exit(EXIT_FAILURE);
		}
	}

	std::string rid("reader_" + std::to_string(1));
	Reader reader_fwd(rid, opts.is_gz);
	vreader.push_back(reader_fwd);
	if (is_paired)
	{
		rid = "reader_" + std::to_string(2);
		Reader reader_rev(rid, opts.is_gz);
		vreader.push_back(reader_rev);
	}

	ss.str("");
	ss << STAMP << "thread: " << std::this_thread::get_id() << " started" << std::endl;
	std::cout << ss.str();
	auto t = std::chrono::high_resolution_clock::now();

	Read read;
	bool done_fwd = false;
	bool done_rev = false;
	uint8_t idx_fwd_reads = 0;
	uint8_t idx_rev_reads = 1;
	// loop calling Readers
	for (; !reader_fwd.is_done || (is_paired && !vreader[idx_rev_reads].is_done);)
	{
		// first push FWD read
		if (!reader_fwd.is_done)
		{
			read = reader_fwd.nextread(ifs_fwd, opts.readfiles[idx_fwd_reads], opts);

			if (!read.isEmpty)
			{
				read.init(opts);
				read.load_db(kvdb); // get matches from Key-value database
				//unmarshallJson(kvdb); // get matches from Key-value database
				readQueue.push(read);
			}
		}
		// second push REV read (if paired)
		if (is_paired && !vreader[1].is_done)
		{
			read = vreader[1].nextread(ifs_rev, opts.readfiles[idx_rev_reads], opts);

			if (!read.isEmpty)
			{
				read.init(opts);
				read.load_db(kvdb); // get matches from Key-value database
				readQueue.push(read);
			}
		}
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
	readQueue.decrPushers(); // signal the reader done adding
	readQueue.notify(); // notify processor that might be waiting to pop

	ss.str("");
	ss << STAMP << "thread: " << std::this_thread::get_id() << " done. Elapsed time: "
		<< std::setprecision(2) << std::fixed << elapsed.count() << " sec Reads added: " << read.id + 1
		<< " readQueue.size: " << readQueue.size() << std::endl;
	std::cout << ss.str();

} // ~ReadControl::run

// ~read_control.cpp
