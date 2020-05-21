#pragma once
/**
 * FILE: reader.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Reads 'Reads file' and pushes the records to a shared queue for further processing
 */

#include <string>
#include <fstream> // std::ifstream

#include "readsqueue.hpp"
//#include "concurrentqueue.h"
#include "kvdb.hpp"
#include "options.hpp"
#include "gzip.hpp"

 // forward
class Read;

struct Readstate {
	Readstate(bool is_gz)
		: is_done(false), isFastq(false), isFasta(false), read_count(0), line_count(0), 
		last_count(0), last_stat(0) {}
	bool is_done;
	bool isFastq; // file is FASTQ
	bool isFasta; // file is FASTA
	unsigned read_count; // count of reads in the file
	unsigned line_count; // count of non-empty lines in the reads file
	unsigned last_count; // count of lines in a single read
	int last_stat;
	std::string last_header; // header line last read
};

/* 
 * reads Reads file and, generates Read objects
 */
class Reader {
public:
	Reader(ReadsQueue& readQueue, std::vector<std::string>& readfiles, bool is_gz);

	void operator()() { run(); }
	void run();

	std::string nextread(std::vector<std::ifstream>& fsl, std::vector<Gzip>& gzips); // new line separated read data. FA - 2 lines, FQ - 4 lines
	bool nextread(const std::string readsfile, std::string &seq);
	void reset();
	static bool hasnext(std::ifstream& ifs);
	static bool loadReadByIdx(Read& read);
	static bool loadReadById(Read& read);

public:
	bool is_done; // flags end of all read streams
	unsigned count_all; // count of reads in all streams

private:
	bool is_gzipped;
	std::size_t next_idx; // index of the reads file to be read next
	std::vector<Readstate> states; // 1st file - FWD, 2dn - REV
	std::vector<std::string>& readfiles;
	ReadsQueue& readQueue;
};

// ~reader.hpp