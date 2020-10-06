#pragma once
/**
 * FILE: Readsfile.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Encapsulates operations on the 'Reads' files both flat and archived
 */

#include <string>
#include <fstream> // std::ifstream

#include "kvdb.hpp"
#include "options.hpp"
#include "izlib.hpp"

 // forward
class Read;

struct Readstate {
	Readstate()	: is_done(false), isFastq(false), isFasta(false), 
		read_count(0), line_count(0), last_count(0), max_reads(0), last_stat(0) {}
	bool is_done;
	bool isFastq; // file is FASTQ
	bool isFasta; // file is FASTA
	unsigned read_count; // count of reads in the file
	unsigned line_count; // count of non-empty lines in the reads file
	unsigned last_count; // count of lines in a single read
	unsigned max_reads;  // max reads expected to be processed
	int last_stat;
	std::string last_header; // header line last read
};

/* 
 * Database like interface for accessing reads
 */
class Readfeed {
public:
	Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, bool is_gz);

	//void operator()() { run(); }
	void run();

	bool next(int fs_idx, std::string& seq); // \n separated read data. FA - 2 lines, FQ - 4 lines
	bool next(std::vector<std::ifstream>& fstreams, std::vector<Izlib>& vzlib, std::vector<Readstate>& vstate, int& inext, std::string& seq);
	std::string nextfwd(std::ifstream& ifs);
	std::string nextrev(std::ifstream& ifs);
	void reset();
	bool split(const unsigned num_parts, const unsigned num_reads, const std::string& outdir);
	static bool is_split_done(const unsigned num_parts, const unsigned num_reads, const std::string dbdir, const std::vector<std::string>& readfiles);
	static bool hasnext(std::ifstream& ifs);
	static bool loadReadByIdx(Read& read);
	static bool loadReadById(Read& read);

public:
	FEED_TYPE type;
	bool is_done; // flags end of all read streams
	bool is_ready; // flags the read feed is ready i.e. no need to run split
	unsigned count_all; // count of reads in all streams

private:
	bool is_gzipped;
	bool is_two_files; // flags two read files are processed (otherwise single file)
	bool is_next_fwd; // flags the next file to be read is FWD (otherwise REV) TODO: remove
	std::vector<std::string>& readfiles;

	// input processing
	std::vector<std::ifstream> ifsv; // [fwd_0, rev_0, fwd_1, rev_1, ... fwd_n-1, rev_n-1]
	std::vector<Izlib> vzlib_in;
	std::vector<Readstate> vstate_in;
	std::vector<bool> vnext_fwd_in; // size = num processors

	// output processing
	std::vector<std::ofstream> ofsv;
	std::vector<Izlib> vzlib_out;
	std::vector<Readstate> vstate_out;
	std::vector<bool> vnext_fwd_out; // size = num processors

	// TODO: remove
	Readstate state_fwd;
	Readstate state_rev;
	Izlib izlib_fwd;
	Izlib izlib_rev;
};

// ~readfeed.hpp