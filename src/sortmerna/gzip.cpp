/* 
 * FILE: gzip.cpp
 * Created: Feb 22, 2018 Thu
 * @copyright 2016-19 Clarity Genomics BVBA
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cassert>
#include <algorithm>

#include "gzip.hpp"


Gzip::Gzip(bool gzipped) 
	: 
	gzipped(gzipped), 
	line_start(0), 
	pstrm(0) 
{ 
	if (gzipped) 
		init(); 
}

/*
 * Called from constructor
 */
void Gzip::init()
{
	static z_stream strm; // TODO: better way i.e. no 'static' ?
	pstrm = &strm;

	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = inflateInit2(&strm, 47);
	if (ret != Z_OK) {
		std::cerr << "Reader::initZstream failed. Error: " << ret << std::endl;
		exit(EXIT_FAILURE);;
	}

	strm.avail_out = 0;

	z_in.resize(IN_SIZE);
	z_out.resize(OUT_SIZE);
	std::fill(z_in.begin(), z_in.end(), 0); // fill IN buffer with 0s
	std::fill(z_out.begin(), z_out.end(), 0); // fill OUT buffer wiht 0s
} // ~Gzip::init

/* 
 * return values: RL_OK (0) | RL_END (1)  | RL_ERR (-1)
 *
 * TODO: Make sure the stream is OK before calling this function.
 *       std::getline doesn't return error if the stream is not 
 *       readable/closed. It returns the same input it was passed.
 */
int Gzip::getline(std::ifstream & ifs, std::string & line)
{
	char* line_end = 0;
	int ret = RL_OK;

	line.clear();

	if (gzipped)
	{
		bool line_ready = false;
		for (; !line_ready; )
		{
			if (!line_start || !line_start[0])
			{
				ret = Gzip::inflatez(ifs); // inflate

				if (ret == Z_STREAM_END && pstrm->avail_out == OUT_SIZE) 
					return RL_END;

				if (ret < 0) 
					return RL_ERR;

				if (!line_start)
					line_start = (char*)z_out.data();
			}

			//line_end = strstr(line_start, "\n"); // returns 0 if '\n' not found
			line_end = std::find(line_start, (char*)&z_out[0] + OUT_SIZE - pstrm->avail_out - 1, 10); // '\n'
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1, [l = std::locale{}](auto ch) { return ch == 10; });
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1, [l = std::locale{}](auto ch) { return std::isspace(ch, l); });
			if (line_end && line_end[0] == 10)
			{
				std::copy(line_start, line_end, std::back_inserter(line));

				if (line_end < (char*)&z_out[0] + OUT_SIZE - pstrm->avail_out - 1) // check there is data after line_end
					line_start = line_end + 1; // skip '\n'
				else
				{
					line_start = 0; // no more data in OUT buffer - flag to inflate more
					pstrm->avail_out = 0; // mark OUT buffer as Full to reflush from the beginning [bug 61]
				}

				line_ready = true; // DEBUG: 
				
				//if (line == "@SRR1635864.196 196 length=101")
				//	std::cout << "HERE";
			}
			else
			{
				line_end = (char*)&z_out[0] + OUT_SIZE - pstrm->avail_out; // end of data in out buffer
				std::copy(line_start, line_end, std::back_inserter(line));
				line_start = (pstrm->avail_out == 0) ? 0 : line_end;
				line_end = 0;
				line_ready = false;
			}
		} // ~for !line_ready
	}
	else // non-compressed file
	{
		if (ifs.eof()) return RL_END;

		std::getline(ifs, line);
		//if (ifs.fail()) return RL_ERR;
	}

	return RL_OK;
} // ~Gzip::getline

/*
 * Called from getline
 */
int Gzip::inflatez(std::ifstream & ifs)
{
	int ret;
	std::stringstream ss;

	for (;;)
	{
		if (pstrm->avail_in == 0 && !ifs.eof()) // in buffer empty
		{
			std::fill(z_in.begin(), z_in.end(), 0); // reset buffer to 0
			ifs.read((char*)z_in.data(), IN_SIZE);
			if (!ifs.eof() && ifs.fail())
			{
				(void)inflateEnd(pstrm);
				return Z_ERRNO;
			}

			pstrm->avail_in = ifs.gcount();
			pstrm->next_in = z_in.data();
		}

		if (pstrm->avail_in == 0 && ifs.eof())
		{
			if (pstrm->avail_out < OUT_SIZE)
				pstrm->avail_out = OUT_SIZE;

			ret = inflateEnd(pstrm); // free up the resources

			if (ret != Z_STREAM_END)
			{
				std::cout << STAMP << "WARNING: inflateEnd status is " << ret << std::endl;
			}

			return Z_STREAM_END;
		}

		if (pstrm->avail_out == 0) // out buffer is full - reset
		{
			std::fill(z_out.begin(), z_out.end(), 0); // reset buffer to 0
			pstrm->avail_out = OUT_SIZE;
			pstrm->next_out = z_out.data();
		}

		ret = inflate(pstrm, Z_NO_FLUSH); //  Z_NO_FLUSH Z_SYNC_FLUSH Z_BLOCK
		assert(ret != Z_STREAM_ERROR);
		switch (ret)
		{
		case Z_NEED_DICT:
			ret = Z_DATA_ERROR; /* and fall through */
		case Z_DATA_ERROR:
		case Z_MEM_ERROR:
			(void)inflateEnd(pstrm);
			return ret;
		case Z_STREAM_END:
			break;
		}

		// OUT buffer holds Inflated data
		// IN buffer holds Compressed data
		// avail_out == 0 means OUT buffer is Full i.e. no space left
		// avail_in  == 0 means  IN buffer is Empty
		// second condition checks if there is still data left in OUT buffer when IN buffer is empty
		if ( pstrm->avail_out == 0 || ( pstrm->avail_out < OUT_SIZE && pstrm->avail_in == 0 ) ) 
			break;
	} // for(;;)

	return ret;// == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
} // ~Gzip::inflatez