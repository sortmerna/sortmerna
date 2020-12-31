#pragma once
/* 
 * FILE: izlib.hpp
 * Created: Feb 22, 2018 Thu
 */

#include <vector>

#include "zlib.h"

#define SIZE_32 32768U /* buffer size 32M */
#define SIZE_16 16384U /* buffer size 16M */
#define RL_OK    0
#define RL_END   1
#define RL_ERR  -1

class Izlib
{
public:
	Izlib(bool is_compress=false, bool is_init=true);

	void init(bool is_compress = false);
	int reset_deflate(); // clean up z_stream
	int finish_deflate(std::ostream& ofs);
	int reset_inflate();
	int getline(std::ifstream & ifs, std::string & line);
	/*
    * deflate passed string and append it to the file stream. Finish processing when the string has 0 size
    * prototype: https://github.com/madler/zlib/blob/master/examples/zpipe.c:def
    * 
    * @param   readstr  a Read as string to be compressed. String of 0 size indicates the end of processing.
    * @param   ofs      compressed output file stream
    * @param   is_last  flags the last string passed -> Finish compressing
    * @return           execution status
    */
	int defstr(std::string& readstr, std::ostream& ofs, bool is_last=false);

private:
	char* line_start; // pointer to the start of a line within the 'z_out' buffer
	z_stream strm; // stream control structure. Holds stream in/out buffers (byte arrays), sizes, positions etc.
	size_t buf_in_size;
	size_t buf_out_size;
	unsigned z_in_num; // number of reads accumulated in IN buffer. For debugging.
	std::vector<unsigned char> z_in; // IN buffer
	std::vector<unsigned char> z_out; // OUT buffer

private:
	int inflatez(std::ifstream & ifs); // 'z' in the name to distinguish from zlib.inflate
};