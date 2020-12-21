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
	int reset_inflate();
	int getline(std::ifstream & ifs, std::string & line);
	int defstr(std::string& readstr, std::ostream& ofs, bool is_last=false);

private:
	char* line_start; // pointer to the start of a line within the 'z_out' buffer
	z_stream strm; // stream control structure. Holds stream in/out buffers (byte arrays), sizes, positions etc.
	size_t buf_in_size;
	size_t buf_out_size;
	std::vector<unsigned char> z_in; // IN buffer
	std::vector<unsigned char> z_out; // OUT buffer

private:
	int inflatez(std::ifstream & ifs); // 'z' in the name to distinguish from zlib.inflate
};