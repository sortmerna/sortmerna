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


Izlib::Izlib(bool gzipped, bool is_compress, bool is_init)
	: 
	gzipped(gzipped), 
	line_start(0),
	strm(),
	buf_in_size(0),
	buf_out_size(0)
{ 
	if (is_init && gzipped) 
		init(is_compress); 
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
void Izlib::init(bool is_compress)
{
	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	int ret = is_compress ? deflateInit(&strm, Z_DEFAULT_COMPRESSION) : inflateInit2(&strm, 47);
	if (ret != Z_OK) {
		ERR("Izlib::init failed. Error: " , ret);
		exit(EXIT_FAILURE);;
	}

	strm.avail_out = 0;

	// compress: in size > out size, inflate: in size < out size
	buf_in_size = is_compress ? SIZE_32 : SIZE_16;
	buf_out_size = is_compress ? SIZE_16 : SIZE_32;
	z_in.resize(buf_in_size);
	z_out.resize(buf_out_size);

	std::fill(z_in.begin(), z_in.end(), 0); // fill IN buffer with 0s
	std::fill(z_out.begin(), z_out.end(), 0); // fill OUT buffer with 0s
} // ~Izlib::init

void Izlib::clean() {
	(void)deflateEnd(&strm);
}

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

				if (ret == Z_STREAM_END && strm.avail_out == buf_out_size)
					return RL_END;

				if (ret < 0) 
					return RL_ERR;

				if (!line_start)
					line_start = (char*)z_out.data();
			}

			//line_end = strstr(line_start, "\n"); // returns 0 if '\n' not found
			line_end = std::find(line_start, (char*)&z_out[0] + buf_out_size - strm.avail_out - 1, 10); // '\n'
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1, [l = std::locale{}](auto ch) { return ch == 10; });
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_SIZE - strm.avail_out - 1, [l = std::locale{}](auto ch) { return std::isspace(ch, l); });
			if (line_end && line_end[0] == 10)
			{
				std::copy(line_start, line_end, std::back_inserter(line));

				if (line_end < (char*)&z_out[0] + buf_out_size - strm.avail_out - 1) // check there is data after line_end
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
				line_end = (char*)&z_out[0] + buf_out_size - strm.avail_out; // end of data in out buffer
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
			ifs.read((char*)z_in.data(), buf_in_size); // get data from reads file into IN buffer 
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
			if (strm.avail_out < buf_out_size)
				strm.avail_out = buf_out_size;

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
			strm.avail_out = buf_out_size;
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
		if ( strm.avail_out == 0 || ( strm.avail_out < buf_out_size && strm.avail_in == 0 ) )
			break;
	} // for(;;)

	return ret;// == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
} // ~Izlib::inflatez

/*
  deflate passed string and append it to the file stream. Finish processing when the string has 0 size
  prototype: https://github.com/madler/zlib/blob/master/examples/zpipe.c:def

  @param   readstr  a Read as string to be compressed. String of 0 size indicates the end of processing.
  @param   ofs      compressed output file stream
  @param   is_last  flags the last string passed -> Finish compressing
  @return           execution status
*/
int Izlib::defstr(std::string& readstr, std::ofstream& ofs, bool is_last)
{
	std::stringstream ss(readstr);
	int ret = Z_OK;
	int flush = Z_NO_FLUSH; // zlib:deflate parameter
	bool is_eos = false; // end of readstr reached
	bool is_deflate = false;
	unsigned pending_bytes = 0;
	int pending_bits = 0;

	// this loop is for repeated reading 'readstr' if
	// it doesn't fit into the IN buffer in one read (very large reads)
	for (; !(is_eos || flush == Z_FINISH || ret == Z_ERRNO);) {
		// add data to IN buffer. Fill up the whole buffer before deflating
		ss.read(reinterpret_cast<char*>(&z_in[0] + strm.avail_in), buf_in_size - strm.avail_in);
		strm.avail_in += ss.gcount();
		if (!ss.eof() && ss.fail()) {
			(void)deflateEnd(&strm);
			ret = Z_ERRNO;
			break;
		}
		
		// deflate or keep accumulating IN?
		if (is_deflate) {
			if (strm.avail_in < buf_in_size && readstr.size() > 0) {
				is_deflate = false;
				break; // keep accumulating IN
			}
		}
		else if (strm.avail_in == buf_in_size || is_last) {
			is_deflate = true; // IN is full - start deflating
		}
		else {
			break; // keep accumulating IN
		}

		strm.next_in = z_in.data();

		flush = is_last ? Z_FINISH : Z_NO_FLUSH; // finish if readstr.size is 0

		// run deflate() until OUT is full i.e. no free space in OUT buffer
		// finish compression if all of source has been read in
		for (;;) {
			if (strm.avail_out == 0) {
				strm.avail_out = buf_out_size;
				strm.next_out = z_out.data();
			}
			ret = deflate(&strm, flush); // runs until OUT is full or IN is empty
			assert(ret != Z_STREAM_ERROR);
			// check accumulated output
			//ret = deflatePending(&strm, &pending_bytes, &pending_bits);
			//assert(ret != Z_STREAM_ERROR);
			// append to the output file (std::ios_base::app)
			ofs.write(reinterpret_cast<char*>(z_out.data()), buf_out_size - strm.avail_out);
			if (ofs.fail()) {
				(void)deflateEnd(&strm);
				ret = Z_ERRNO;
				break;
			}
			if (strm.avail_out > 0 || (strm.avail_out == 0 && is_last))
				break;
		} // ~for

		assert(strm.avail_in == 0); // all input was used
		is_eos = ss.eof();
	} // ~for

	if (is_deflate && !is_last) {
		assert(ret == Z_STREAM_END); // stream will be complete
	}

	if (flush == Z_FINISH) {
		(void)deflateEnd(&strm);
	}
	return ret;
} // ~Izlib::deflatez