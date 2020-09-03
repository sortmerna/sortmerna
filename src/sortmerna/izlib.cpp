/* 
 * FILE: gzip.cpp
 * Created: Feb 22, 2018 Thu
 * @copyright 2016-20 Clarity Genomics BVBA
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <cassert>
#include <algorithm>

#include "izlib.hpp"
#include "common.hpp"


Izlib::Izlib(bool gzipped)
	: 
	gzipped(gzipped), 
	line_start(0)
{ 
	if (gzipped) 
		init(); 
}


//Izlib::~Izlib() {
//	line_start = 0;
//	strm.zalloc = Z_NULL;
//	strm.zfree = Z_NULL;
//	strm.opaque = Z_NULL;
//	strm.avail_in = Z_NULL;
//	strm.next_in = Z_NULL;
//	strm.avail_out = Z_NULL;
//}

/*
 * Called from constructor
 */
void Izlib::init()
{
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = inflateInit2(&strm, 47);
	if (ret != Z_OK) {
		ERR("Reader::initZstream failed. Error: " , ret);
		exit(EXIT_FAILURE);;
	}

	strm.avail_out = 0;

	z_in.resize(IN_SIZE);
	z_out.resize(OUT_SIZE);
	std::fill(z_in.begin(), z_in.end(), 0); // fill IN buffer with 0s
	std::fill(z_out.begin(), z_out.end(), 0); // fill OUT buffer with 0s
} // ~Gzip::init

/* 
  get a line from the compressed stream
 
  TODO: Make sure the stream is OK before calling this function.
        std::getline doesn't return error if the stream is not 
        readable/closed. It returns the same input it was passed.

  return values: RL_OK (0) | RL_END (1)  | RL_ERR (-1)
 */
int Izlib::getline(std::ifstream& ifs, std::string& line)
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
				ret = Izlib::inflatez(ifs); // inflate

				if (ret == Z_STREAM_END && strm.avail_out == OUT_SIZE) 
					return RL_END;

				if (ret < 0) 
					return RL_ERR;

				if (!line_start)
					line_start = (char*)z_out.data();
			}

			//line_end = strstr(line_start, "\n"); // returns 0 if '\n' not found
			line_end = std::find(line_start, (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1, 10); // '\n'
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1, [l = std::locale{}](auto ch) { return ch == 10; });
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1, [l = std::locale{}](auto ch) { return std::isspace(ch, l); });
			if (line_end && line_end[0] == 10)
			{
				std::copy(line_start, line_end, std::back_inserter(line));

				if (line_end < (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1) // check there is data after line_end
					line_start = line_end + 1; // skip '\n'
				else
				{
					line_start = 0; // no more data in OUT buffer - flag to inflate more
					strm.avail_out = 0; // mark OUT buffer as Full to reflush from the beginning [bug 61]
				}

				line_ready = true; // DEBUG: 
				
				//if (line == "@SRR1635864.196 196 length=101")
				//	std::cout << "HERE";
			}
			else
			{
				line_end = (char*)&z_out[0] + OUT_SIZE - strm.avail_out; // end of data in out buffer
				std::copy(line_start, line_end, std::back_inserter(line));
				line_start = (strm.avail_out == 0) ? 0 : line_end;
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
} // ~Izlib::getline

/*
 * Called from getline
 */
int Izlib::inflatez(std::ifstream & ifs)
{
	int ret;

	for (;;)
	{
		if (strm.avail_in == 0 && !ifs.eof()) // in buffer empty
		{
			std::fill(z_in.begin(), z_in.end(), 0); // reset IN buffer to 0
			ifs.read((char*)z_in.data(), IN_SIZE); // get data from reads file into IN buffer 
			if (!ifs.eof() && ifs.fail()) // not end of reads file And read fail -> round up and return error
			{
				(void)inflateEnd(&strm);
				return Z_ERRNO;
			}

			strm.avail_in = ifs.gcount();
			strm.next_in = z_in.data();
		}

		if (strm.avail_in == 0 && ifs.eof())
		{
			if (strm.avail_out < OUT_SIZE)
				strm.avail_out = OUT_SIZE;

			ret = inflateEnd(&strm); // free up the resources

			if (ret != Z_STREAM_END)
			{
				INFO("inflateEnd status is ", ret);
			}

			return Z_STREAM_END;
		}

		if (strm.avail_out == 0) // out buffer is full - reset
		{
			std::fill(z_out.begin(), z_out.end(), 0); // reset buffer to 0
			strm.avail_out = OUT_SIZE;
			strm.next_out = z_out.data();
		}

		ret = inflate(&strm, Z_NO_FLUSH); //  Z_NO_FLUSH Z_SYNC_FLUSH Z_BLOCK
		assert(ret != Z_STREAM_ERROR);
		switch (ret)
		{
		case Z_NEED_DICT:
			ret = Z_DATA_ERROR; /* and fall through */
		case Z_DATA_ERROR:
		case Z_MEM_ERROR:
			(void)inflateEnd(&strm);
			return ret;
		case Z_STREAM_END:
			break;
		}

		// OUT buffer holds Inflated data
		// IN buffer holds Compressed data
		// avail_out == 0 means OUT buffer is Full i.e. no space left
		// avail_in  == 0 means  IN buffer is Empty
		// second condition checks if there is still data left in OUT buffer when IN buffer is empty
		if ( strm.avail_out == 0 || ( strm.avail_out < OUT_SIZE && strm.avail_in == 0 ) ) 
			break;
	} // for(;;)

	return ret;// == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
} // ~Izlib::inflatez

/*
  deflate passed string and append it to the output file
  prototype: https://github.com/madler/zlib/blob/master/examples/zpipe.c:def

  @param   readstr  read as string
  @param   ofs      compressed output file stream
  @return           execution status
*/
int Izlib::deflatez(std::string& readstr, std::ofstream& ofs)
{
	std::stringstream readss(readstr);
	int ret = Z_OK;
	int flush = Z_NO_FLUSH; // zlib:deflate parameter

	for (; flush != Z_FINISH && ret != Z_ERRNO;) {
		readss.read(reinterpret_cast<char*>(&z_in[0]), IN_SIZE);
		strm_def.avail_in = readss.gcount();
		if (!readss.eof() && readss.fail()) {
			(void)deflateEnd(&strm_def);
			ret = Z_ERRNO;
			break;
		}
		flush = readss.eof() ? Z_FINISH : Z_NO_FLUSH;
		strm_def.next_in = z_in.data();

		// run deflate() on input until output buffer not full,
		// finish compression if all of source has been read in
		for (; strm_def.avail_out == 0;) {
			strm_def.avail_out = OUT_SIZE;
			strm_def.next_out = z_out.data();
			ret = deflate(&strm_def, flush);
			assert(ret != Z_STREAM_ERROR);
			// append to the output file (std::ios_base::app)
			ofs.write(reinterpret_cast<char*>(z_out.data()), OUT_SIZE - strm_def.avail_out);
			if (ofs.fail()) {
				(void)deflateEnd(&strm_def);
				ret = Z_ERRNO;
				break;
			}
			// done when last data in file processed
		} // ~for

		assert(strm_def.avail_in == 0); // all input will be used
	} // ~for

	assert(ret == Z_STREAM_END); // stream will be complete

	(void)deflateEnd(&strm_def); // clean up
	return Z_OK;
} // ~Izlib::deflatez