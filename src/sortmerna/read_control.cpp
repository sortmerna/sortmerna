/**
 * FILE: read_control.cpp 
 * Created: Mar 12, 2019 Tue
 *
 * @copyright 2016-20 Clarity Genomics BVBA
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

	size_t read_cnt = 0;
	std::size_t num_aligned = 0; // count of aligned reads (passing E-value)
	bool done_fwd = false;
	bool done_rev = false;
	uint8_t IDX_FWD_READS = 0;
	uint8_t IDX_REV_READS = 1;

	bool is_two_reads = opts.readfiles.size() == 2; // i.e. 2 read files are supplied

	// init FWD Reader
	Reader reader_fwd("reader_fwd", opts.is_gz);
	auto fwd_file = opts.readfiles[IDX_FWD_READS];
	std::ifstream ifs_fwd(fwd_file, std::ios_base::in | std::ios_base::binary);

	if (!ifs_fwd.is_open()) 
	{
		ERR("failed to open file: [" + fwd_file + "]");
		exit(EXIT_FAILURE);
	}

	// init REV Reader
	std::ifstream ifs_rev;
	Reader reader_rev("reader_rev", opts.is_gz);
	if (is_two_reads)
	{
		auto rev_file = opts.readfiles[IDX_REV_READS];
		ifs_rev.open(rev_file, std::ios_base::in | std::ios_base::binary);

		if (!ifs_rev.is_open()) {
			ERR("failed to open file: [" + rev_file + "]");
			exit(EXIT_FAILURE);
		}
	}

	ss.str("");
	ss << STAMP << "thread: " << std::this_thread::get_id() << " started" << std::endl;
	std::cout << ss.str();
	auto t = std::chrono::high_resolution_clock::now();

	// loop calling Readers
	for (; !reader_fwd.is_done || (is_two_reads && !reader_rev.is_done);)
	{
		// first push FWD read
		if (!reader_fwd.is_done)
		{
			Read read = reader_fwd.nextread(ifs_fwd, IDX_FWD_READS, opts);

			if (!read.isEmpty)
			{
				read.init(opts);
				read.load_db(kvdb); // get matches from Key-value database
				//unmarshallJson(kvdb); // get matches from Key-value database
				++read_cnt; // save because push(read) uses move(read)
				if (read.hit) ++num_aligned;
				readQueue.push(read);
			}
		}
		// second push REV read (if paired)
		if (is_two_reads && !reader_rev.is_done)
		{
			Read read = reader_rev.nextread(ifs_rev, IDX_REV_READS, opts);

			if (!read.isEmpty)
			{
				read.init(opts);
				read.load_db(kvdb); // get matches from Key-value database
				++read_cnt;
				if (read.hit) ++num_aligned;
				readQueue.push(read);
			}
		}
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
	readQueue.decrPushers(); // signal the reader done adding
	readQueue.notify(); // notify processor that might be waiting to pop

	ss.str("");
	ss << STAMP << "thread: " << std::this_thread::get_id() << " done. Elapsed time: "
		<< std::setprecision(2) << std::fixed << elapsed.count() << " sec Reads added: " << read_cnt
		<< " Num aligned reads (passing E-value): " << num_aligned
		<< " readQueue.size: " << readQueue.size() << std::endl;
	std::cout << ss.str();

} // ~ReadControl::run

// ~read_control.cpp
