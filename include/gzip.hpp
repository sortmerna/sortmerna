#pragma once
/* 
 * FILE: gzip.hpp
 * Created: Feb 22, 2018 Thu
 */

#include <vector>

#include "zlib.h"

#include "options.hpp"

#define OUT_CHUNK 32768U /* out buffer size */
#define CHUNK 16384      /* file input buffer size */
#define RL_OK 0
#define RL_END 1
#define RL_ERR -1

class Gzip
{
public:
	Gzip(Runopts & opts) : opts(opts), line_start(0), pstrm(0) { if (opts.have_reads_gz) init(); }

	int getline(std::ifstream & ifs, std::string & line);

private:
	Runopts & opts;
	// zlib related
	char* line_start;
	z_stream * pstrm;
	std::vector<unsigned char> z_in;
	std::vector<unsigned char> z_out;

private:
	void init();
	int inf(std::ifstream & ifs);
};