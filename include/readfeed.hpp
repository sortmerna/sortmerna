#pragma once
/**
 * FILE: Readsfile.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Encapsulates operations on the 'Reads' files both flat and archived
 */

#include <string>
#include <fstream> // std::ifstream
#include <filesystem>

#include "kvdb.hpp"
#include "options.hpp"
#include "izlib.hpp"

 // forward
class Read;

struct Readstate {
	Readstate()	: is_done(false), read_count(0), line_count(0), last_count(0), last_stat(0) {}
	void reset() { is_done = false; read_count = 0; line_count = 0; 
		last_count = 0; last_stat = 0; last_header.clear(); }
	bool is_done; // flags EOF reached
	unsigned read_count; // count of reads in the file
	unsigned line_count; // count of non-empty lines in the reads file
	unsigned last_count; // count of lines in a single read
	int last_stat;
	std::string last_header; // header line last read
};

struct Readfile {
	Readfile() : isFastq(false), isFasta(false), isZip(false), size(0), numreads(0) {}
	bool isFastq; // file is FASTQ
	bool isFasta; // file is FASTA
	bool isZip;   // true (compressed) | false (flat)
	unsigned numreads;  // max reads expected to be processed
	std::filesystem::path path;
	std::streampos size;
};

/* 
 * Database like interface for accessing reads
 */
class Readfeed {
public:
	Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, std::filesystem::path& basedir);
	Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, const unsigned num_parts, std::filesystem::path& basedir);

	void run();
	bool next(int inext, std::string& readstr);
	void reset();
	void rewind();
	void rewind_in();
	void init(std::vector<std::string>& readfiles);
	void init_split_files();
	void init_reading();
	void init_vzlib_in();
	//void init_vstate_in();
	bool split();
	bool define_format();
	void count_reads();
	void write_descriptor();
	bool is_split_ready();
	static bool is_split_done(const unsigned num_parts, const unsigned num_reads, const std::string dbdir, const std::vector<std::string>& readfiles);
	static bool hasnext(std::ifstream& ifs);
	static bool loadReadByIdx(Read& read);
	static bool loadReadById(Read& read);

private:
	bool next(int inext, std::string& readstr, bool is_orig, std::vector<Readfile>& files);
	//bool next(int inext, std::string& readstr, unsigned& readlen, bool is_orig = false); // \n separated read data. FA - 2 lines, FQ - 4 lines
	bool next(int inext, std::string& readstr, unsigned& readlen, bool is_orig, std::vector<Readfile>& files);

public:
	FEED_TYPE type;
	BIO_FORMAT biof;
	ZIP_FORMAT zipf;
	bool is_done; // flags end of all read streams
	bool is_ready; // flags the read feed is ready i.e. no need to run split
	bool is_format_defined; // flags the file format is defined i.e. 'define_format' was success
	bool is_two_files; // flags two read files are processed (otherwise single file)
	unsigned num_orig_files; // number of original reads files
	unsigned num_splits;
	unsigned num_reads_tot; // count of reads in all streams
	unsigned length_all; // length of all reads from all files
	unsigned min_read_len;
	unsigned max_read_len;
	std::filesystem::path& basedir; // split files root directory opts.readb
	std::vector<Readfile> orig_files;
private:
	std::vector<Readfile> split_files;

	// input processing
	std::vector<std::ifstream> ifsv; // [fwd_0, rev_0, fwd_1, rev_1, ... fwd_n-1, rev_n-1]
	std::vector<Izlib> vzlib_in;
	std::vector<Readstate> vstate_in;

	// output processing
	std::vector<std::ofstream> ofsv;
	std::vector<Izlib> vzlib_out;
	std::vector<Readstate> vstate_out;
};

// ~readfeed.hpp