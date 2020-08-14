#pragma once
/* 
 * FILE: gzip.hpp
 * Created: Feb 22, 2018 Thu
 */

#include <vector>

#include "zlib.h"

#include "options.hpp"

#define OUT_SIZE 32768U /* out buffer size */
#define  IN_SIZE 16384  /* file input buffer size */
#define    RL_OK 0
#define   RL_END 1
#define   RL_ERR -1

class Izlib
{
public:
	Izlib(bool gzipped);
	//~Izlib();

	void init();
	int getline(std::ifstream & ifs, std::string & line);

private:
	bool gzipped;
	// zlib related
	char* line_start; // pointer to the start of a line within the 'z_out' buffer
	z_stream strm; // stream control structure. Holds stream in/out buffers (byte arrays), sizes, positions etc.
	std::vector<unsigned char> z_in; // IN buffer for compressed data
	std::vector<unsigned char> z_out; // OUT buffer for decompressed data

private:
	int inflatez(std::ifstream & ifs); // 'z' in the name to distinguish from zlib.inflate
};