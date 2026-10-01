/*
@copyright 2016-2026 Clarity Genomics Inc
@copyright 2012-2016 Bonsai Bioinformatics Research Group
@copyright 2014-2016 Knight Lab, Department of Pediatrics, UCSD, La Jolla

@parblock
SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA

This is a free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SortMeRNA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
@endparblock

@contributors Jenya Kopylova   jenya.kopylov@gmail.com
              Laurent Noé      laurent.noe@lifl.fr
              Pierre Pericard  pierre.pericard@lifl.fr
              Daniel McDonald  wasade@gmail.com
              Mikaël Salson    mikael.salson@lifl.fr
              Hélène Touzet    helene.touzet@lifl.fr
              Rob Knight       robknight@ucsd.edu
*/

/* 
 * FILE: gzip.cpp
 * Created: Feb 22, 2018 Thu
 */

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include "izlib.hpp"
#include "common.hpp"

namespace sortmerna {

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
		SMR_THROW("Izlib::init failed. Error: " , ret);
	}
	strm_active = true;

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
	if (!strm_active) return Z_OK;
	strm_active = false;
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
		if (ret == Z_STREAM_ERROR)
			SMR_THROW("Z_STREAM_ERROR in reset_inflate - zlib stream is inconsistent");

		if (ret == Z_STREAM_END) {
			ret = inflateReset(&strm);
			if (ret != Z_OK)
				SMR_THROW("inflateReset failed in reset_inflate: ", ret);
		}
		else if (ret != Z_OK) {
			break; // e.g. Z_DATA_ERROR on trailing garbage: nothing left to drain
		}
	}
	if (strm_active) {
		strm_active = false;
		ret = inflateEnd(&strm); // ret is Z_OK, not Z_STREAM_END
		if (ret != Z_OK && ret != Z_STREAM_END)
			SMR_THROW("inflateEnd failed in reset_inflate: ", ret);
	}
	return ret;
}

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

int Izlib::inflatez(std::ifstream& ifs)
{
	int ret;

	for (;;)
	{
		// IN buffer empty - fill up with more data
		if (strm.avail_in == 0 && !ifs.eof()) // in buffer empty
		{
			std::fill(z_in.begin(), z_in.end(), 0); // reset IN buffer to 0
			ifs.read((char*)z_in.data(), buf_in_size); // add data into IN buffer 
			if (!ifs.eof() && ifs.fail()) // not end of reads file And read fail -> round up and return error
			{
				strm_active = false;
				inflateEnd(&strm);
				return Z_ERRNO;
			}

			strm.avail_in = ifs.gcount();
			strm.next_in = z_in.data();
		}

		// EOF reached - return
		if (strm.avail_in == 0 && ifs.eof())
		{
			if (strm.avail_out < buf_out_size) strm.avail_out = buf_out_size;
			if (strm_active) {
				strm_active = false;
				ret = inflateEnd(&strm); // free up the resources
				if (ret != Z_OK)
					SMR_THROW("inflateEnd failed at EOF in inflatez: ", ret);
			}
			return Z_STREAM_END;
		}

		// OUT buffer is full - reset
		if (strm.avail_out == 0)
		{
			std::fill(z_out.begin(), z_out.end(), 0); // reset buffer to 0
			strm.avail_out = buf_out_size;
			strm.next_out = z_out.data();
		}

		ret = inflate(&strm, Z_NO_FLUSH); //  Z_NO_FLUSH Z_SYNC_FLUSH Z_BLOCK
		if (ret == Z_STREAM_ERROR)
			SMR_THROW("Z_STREAM_ERROR in inflatez - zlib stream is inconsistent");

		switch (ret)
		{
		case Z_NEED_DICT:
			ret = Z_DATA_ERROR; /* and fall through */
		case Z_DATA_ERROR:
		case Z_MEM_ERROR:
			if (strm_active) { strm_active = false; inflateEnd(&strm); }
			return ret;
		case Z_STREAM_END:
			ret = inflateReset(&strm);
			if (ret != Z_OK)
				SMR_THROW("inflateReset failed in inflatez: ", ret);
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

int Izlib::defstr(const std::string& readstr, std::ostream& ofs, bool is_last, const int&& dbg)
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
			strm_active = false;
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
			if (ret == Z_STREAM_ERROR)
				SMR_THROW("Z_STREAM_ERROR in defstr - zlib stream is inconsistent");
			// check accumulated output
			//ret = deflatePending(&strm, &pending_bytes, &pending_bits);
			//assert(ret != Z_STREAM_ERROR);
			// append to the output file (std::ios_base::app)
			ofs.write(reinterpret_cast<char*>(z_out.data()), buf_out_size - strm.avail_out);
			if (ofs.fail()) {
				strm_active = false;
				deflateEnd(&strm);
				ret = Z_ERRNO;
				break;
			}
			if (flush == Z_FINISH && (ret == Z_OK || ret == Z_BUF_ERROR))
				continue;
			if (strm.avail_out > 0 || (strm.avail_out == 0 && is_last))
				break;
		} // ~for

		if (strm.avail_in != 0) // all input must be used
			SMR_THROW("Unexpected remaining input after deflate in defstr");
		if (is_ess) z_in_num = 0; // reset
	} // ~for

	if (flush == Z_FINISH) {
		if (ret != Z_STREAM_END)
			SMR_THROW("Expected Z_STREAM_END after final deflate in defstr, got: ", ret);
		strm_active = false;
		deflateEnd(&strm);
		ofs.flush();
		if (dbg > 1)
			INFO("deflateEnd called");
	}
	return ret;
} // ~Izlib::defstr

/*
* flush the stream and reset deflate
*/
int Izlib::finish_deflate(std::ostream& ofs, const int&& dbg)
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
		if (ret == Z_STREAM_ERROR)
			SMR_THROW("Z_STREAM_ERROR in finish_deflate - zlib stream is inconsistent");
		ofs.write(reinterpret_cast<char*>(z_out.data()), buf_out_size - strm.avail_out);
		if (ofs.fail()) {
			strm_active = false;
			(void)deflateEnd(&strm);
			ret = Z_ERRNO;
			break;
		}
		if (ret == Z_OK || ret == Z_BUF_ERROR)
			continue;
		if (ret == Z_STREAM_END) break;
	} // ~for
	ofs.flush();
	if (dbg > 1)
		INFO("deflateEnd called");
	strm_active = false;
	return deflateEnd(&strm);
}

} // namespace sortmerna