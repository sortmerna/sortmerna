#pragma once

#include <string>

/*
 * reading state of a file
 */
struct Readstate {
	Readstate() : is_done(false), read_count(0), line_count(0), last_count(0), last_stat(0), last_header("") {}
	void reset() {
		is_done = false; read_count = 0; line_count = 0;
		last_count = 0; last_stat = 0; last_header.clear();
	}
	bool is_done; // flags EOF reached
	unsigned read_count; // count of reads in the file
	unsigned line_count; // count of non-empty lines in the reads file
	unsigned last_count; // count of lines in a single read
	int last_stat;
	std::string last_header; // header line last read
};
