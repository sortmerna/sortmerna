#pragma once
/* 
 * FILE: izlib.hpp
 * Created: Feb 22, 2018 Thu
 */

#include <vector>

#include "zlib.h"

#define OUT_SIZE 32768U /* out buffer size */
#define IN_SIZE  16384  /* file input buffer size */
#define RL_OK    0
#define RL_END   1
#define RL_ERR  -1

class Izlib
{
public:
	/*
	  @param gzipped      flags the file is gzipped (true) otherwise flat (false)
	  @param is_compress  flags to compress (true) or inflate (false) the output
	  @param is_init      flags to perform the zlib stream initialization
	*/
	Izlib(bool gzipped=true, bool is_compress=false, bool is_init=true);

	void init(bool is_compress = false);
	void clean(); // clean up z_stream
	int getline(std::ifstream & ifs, std::string & line);
	int deflatez(std::string& readstr, std::ofstream& ofs);

private:
	bool gzipped;
	char* line_start; // pointer to the start of a line within the 'z_out' buffer
	z_stream strm; // stream control structure. Holds stream in/out buffers (byte arrays), sizes, positions etc.
	std::vector<unsigned char> z_in; // IN buffer
	std::vector<unsigned char> z_out; // OUT buffer

private:
	int inflatez(std::ifstream & ifs); // 'z' in the name to distinguish from zlib.inflate
};