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


/*
 * @param is_compress  flags to compress (true) or inflate (false) the output
 * @param is_init      flags to perform the zlib stream initialization
 */
Izlib::Izlib(bool is_compress, bool is_init)
	: 
	line_start(0),
	strm(),
	buf_in_size(0),
	buf_out_size(0),
	z_in_num(0)
{ 
	if (is_init) 
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
	int windowBits = 15;
	int GZIP_ENCODING = 16;
	int memLevel = 8; // default value see zlib.h

	strm.zalloc = Z_NULL;
	strm.zfree = Z_NULL;
	strm.opaque = Z_NULL;
	strm.avail_in = 0;
	strm.next_in = Z_NULL;
	// use deflateInit2 and inflateInit2 to output/input gzip format
	int ret = is_compress 
		? deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED, windowBits | GZIP_ENCODING, memLevel, Z_DEFAULT_STRATEGY)
		: inflateInit2(&strm, 47);
	if (ret != Z_OK) {
		ERR("Izlib::init failed. Error: " , ret);
		exit(EXIT_FAILURE);;
	}

	strm.avail_out = 0;

	line_start = 0;
	// compress: in size > out size, inflate: in size < out size
	buf_in_size = is_compress ? SIZE_32 : SIZE_16;
	buf_out_size = is_compress ? SIZE_16 : SIZE_32;
	z_in.resize(buf_in_size);
	z_out.resize(buf_out_size);

	std::fill(z_in.begin(), z_in.end(), 0); // fill IN buffer with 0s
	std::fill(z_out.begin(), z_out.end(), 0); // fill OUT buffer with 0s
} // ~Izlib::init

int Izlib::reset_deflate() {
	return deflateEnd(&strm);
}

/*
  TODO: 20201013: inflate(Z_FINISH) gives Z_BUF_ERROR (-5) always. inflate(Z_NO_FLUSH) gives Z_OK. Why?
*/
int Izlib::reset_inflate()
{
	int ret = 0;
	for (;strm.avail_in > 0;) {
		std::fill(z_out.begin(), z_out.end(), 0); // reset buffer to 0
		strm.avail_out = buf_out_size;
		strm.next_out = &z_out[0];
		ret = inflate(&strm, Z_NO_FLUSH); // Z_FINISH -> Z_BUF_ERROR, Z_NO_FLUSH -> Z_OK
		assert(ret == Z_OK || ret == Z_STREAM_END);
	}
	ret = inflateEnd(&strm); // ret is Z_OK, not Z_STREAM_END
	assert(ret == Z_OK || ret == Z_STREAM_END);
	return ret;
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
			ifs.read((char*)z_in.data(), buf_in_size); // add data into IN buffer 
			if (!ifs.eof() && ifs.fail()) // not end of reads file And read fail -> round up and return error
			{
				inflateEnd(&strm);
				return Z_ERRNO;
			}

			strm.avail_in = ifs.gcount();
			strm.next_in = z_in.data();
		}

		if (strm.avail_in == 0 && ifs.eof())
		{
			if (strm.avail_out < buf_out_size) strm.avail_out = buf_out_size;
			ret = inflateEnd(&strm); // free up the resources
			assert(ret == Z_OK || ret == Z_STREAM_END); // seems always Z_OK, never Z_STREAM_END - why?
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
			inflateEnd(&strm);
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

int Izlib::defstr(std::string& readstr, std::ostream& ofs, bool is_last)
{
	std::stringstream ss(readstr);
	int ret = Z_OK;
	int flush = Z_NO_FLUSH; // zlib:deflate parameter
	bool is_ess = false; // end of readstr reached
	bool is_deflate = false; // flags to run deflation
	//unsigned pending_bytes = 0;
	//int pending_bits = 0;

	// this loop is for repeated reading 'readstr' if
	// it doesn't fit into the IN buffer in one read
	for (; !(is_ess || flush == Z_FINISH || ret == Z_ERRNO);) {
		// add data to IN buffer. Fill up the whole buffer before deflating
		if (strm.avail_in == 0) strm.next_in = &z_in[0];
		ss.read(reinterpret_cast<char*>(&z_in[0] + strm.avail_in), buf_in_size - strm.avail_in);
		strm.avail_in += ss.gcount();
		if (!ss.eof() && ss.fail()) {
			deflateEnd(&strm);
			ret = Z_ERRNO;
			break;
		}
		is_ess = ss.eof();
		if (is_ess) ++z_in_num;
		
		// if IN not full and not last read -> return for more reads
		if (is_deflate) {
			if (strm.avail_in < buf_in_size && !is_last) {
				is_deflate = false;
				break; // keep accumulating IN
			}
		}
		else if (strm.avail_in == buf_in_size || is_last) {
			is_deflate = true; // IN is full, or not full and last string - deflate
		}
		else {
			break; // keep accumulating IN
		}

		flush = is_last && is_ess ? Z_FINISH : Z_NO_FLUSH; // finish if the last string + end of the last string

		// run deflate() until OUT is full i.e. no free space in OUT buffer
		// finish compression if all of source has been read in
		for (;;) {
			// reset the OUT buffer as all deflated output was written to file at this point
			//std::fill(z_out.begin(), z_out.end(), 0); // reset buffer to 0
			strm.avail_out = buf_out_size;
			strm.next_out = z_out.data();
			// deflate
			ret = deflate(&strm, flush); // runs until OUT is full or IN is empty
			assert(ret != Z_STREAM_ERROR);
			// check accumulated output
			//ret = deflatePending(&strm, &pending_bytes, &pending_bits);
			//assert(ret != Z_STREAM_ERROR);
			// append to the output file (std::ios_base::app)
			ofs.write(reinterpret_cast<char*>(z_out.data()), buf_out_size - strm.avail_out);
			if (ofs.fail()) {
				deflateEnd(&strm);
				ret = Z_ERRNO;
				break;
			}
			if (flush == Z_FINISH && (ret == Z_OK || ret == Z_BUF_ERROR))
				continue;
			if (strm.avail_out > 0 || (strm.avail_out == 0 && is_last))
				break;
		} // ~for

		assert(strm.avail_in == 0); // all input was used
		if (is_ess) z_in_num = 0; // reset
	} // ~for

	if (flush == Z_FINISH) {
		assert(ret == Z_STREAM_END);
		deflateEnd(&strm);
		ofs.flush();
		INFO("deflateEnd called");
	}
	return ret;
} // ~Izlib::defstr

/*
* flush the stream and reset deflate
*/
int Izlib::finish_deflate(std::ostream& ofs)
{
	int ret = Z_OK;
	// run deflate() until OUT is full i.e. no free space in OUT buffer
	// finish compression if all of source has been read in
	for (;;) {
		// reset the OUT buffer as all deflated output was written to file at this point
		//std::fill(z_out.begin(), z_out.end(), 0); // reset buffer to 0
		strm.avail_out = buf_out_size;
		strm.next_out = z_out.data();
		// deflate
		ret = deflate(&strm, Z_FINISH); // runs until OUT is full or IN is empty
		assert(ret != Z_STREAM_ERROR);
		ofs.write(reinterpret_cast<char*>(z_out.data()), buf_out_size - strm.avail_out);
		if (ofs.fail()) {
			(void)deflateEnd(&strm);
			ret = Z_ERRNO;
			break;
		}
		if (ret == Z_OK || ret == Z_BUF_ERROR)
			continue;
		if (ret == Z_STREAM_END) break;
	} // ~for
	ofs.flush();
	INFO("deflateEnd called");
	return deflateEnd(&strm);
}