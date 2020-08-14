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
#include "kvdb.hpp"
#include "options.hpp"
#include "izlib.hpp"

 // forward
class Read;

struct Readstate {
	Readstate()	: is_done(false), isFastq(false), isFasta(false), 
		read_count(0), line_count(0), last_count(0), last_stat(0) {}
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
class Readsfile {
public:
	Readsfile(ReadsQueue& readQueue, std::vector<std::string>& readfiles, bool is_gz);

	void operator()() { run(); }
	void run();

	std::string nextread(std::ifstream& ifs); // new line separated read data. FA - 2 lines, FQ - 4 lines
	std::string nextfwd(std::ifstream& ifs);
	std::string nextrev(std::ifstream& ifs);
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
	bool is_two_files; // flags two read files are processed (otherwise single file)
	bool is_next_fwd; // flags the next file to be read is FWD (otherwise REV)
	std::vector<std::string>& readfiles;
	ReadsQueue& readQueue;
	Readstate state_fwd;
	Readstate state_rev;
	Izlib izlib_fwd;
	Izlib izlib_rev;
};

// ~reader.hpp