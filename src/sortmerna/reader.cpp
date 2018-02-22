/**
 * FILE: reader.cpp
 * Created: Nov 26, 2017 Sun
 * 
 * Processes Reads file, creates Read objects, and pushes them to a queue for further pick−up by Processor
 */

#include <string>
#include <locale> // std::isspace
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <chrono> // std::chrono
#include <ios> // std::ios_base
#include <iomanip> // std::precision
#include <vector>
#include <algorithm> // find, find_if
//#include <cctype> // isspace

#include <zlib.h>

#include "reader.hpp"

void Reader::read()
{
	std::stringstream ss;

	std::ifstream ifs(opts.readsfile, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) {
		ss << "failed to open " << opts.readsfile << std::endl;
		std::cout << ss.str(); ss.str("");
		exit(EXIT_FAILURE);
	}
	else
	{
		std::string line;
		Read read; // an empty read
		unsigned int read_id = 0; // read ID
		bool isFastq = true;
		bool lastRec = false; // lastRec is to make one iteration past the EOF

		if (opts.have_reads_gz)
		{
			initZstream();
		}

		ss << id << " thread: " << std::this_thread::get_id() << " started\n";
		std::cout << ss.str(); ss.str("");
		auto t = std::chrono::high_resolution_clock::now();

		for (int count = 0; ; ) // count lines in a single record
		{
			if (!lastRec) getline(ifs, line); // std::getline(ifs, line)

			if (line.empty() && !lastRec)
			{
				if (ifs.eof()) lastRec = true;
				continue;
			}

			if (lastRec)
			{
				if (!read.isEmpty)
				{
					read.init(opts, kvdb, read_id); // load alignment statistics from DB
					readQueue.push(read);
					++read_id;
				}
				break;
			}

			// left trim space and '>' or '@'
			//line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));
			// right-trim whitespace in place (removes '\r' too)
			line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
			// removes all space
			//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));

			// fastq: 0(header), 1(seq), 2(+), 3(quality)
			// fasta: 0(header), 1(seq)
			if (line[0] == FASTA_HEADER_START || line[0] == FASTQ_HEADER_START)
			{
				if (!read.isEmpty)
				{
					read.init(opts, kvdb, read_id);
					readQueue.push(read);
					++read_id;
					count = 0;
				}

				// start new record
				read.clear();
				isFastq = (line[0] == FASTQ_HEADER_START);
				read.format = isFastq ? Format::FASTQ : Format::FASTA;
				read.header = line;
				read.isEmpty = false;
			} // ~if header line
			else {
				++count;
				if (isFastq && line[0] == '+') continue;
				if (isFastq && count == 3)
				{
					read.quality = line;
					continue;
				}
				read.sequence += line;
			}
			if (ifs.eof()) lastRec = true; // push and break
		} // ~for getline

		std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
		readQueue.mDoneAdding();
		ss << id << " thread: " << std::this_thread::get_id() << " done. Elapsed time: " 
			<< std::setprecision(2) << std::fixed << elapsed.count() << " sec Reads added: " << read_id << std::endl;
		std::cout << ss.str(); ss.str("");
	}
	ifs.close();
} // ~Reader::read

int Reader::getline(std::ifstream & ifs, std::string & line)
{
	char* line_end = 0;
	int ret = RL_OK;

	if (opts.have_reads_gz)
	{
		bool line_ready = false;
		for (; !line_ready; )
		{
			if (!line_start || !line_start[0])
			{
				ret = inf(ifs); // inflate

				if (ret == Z_STREAM_END && pstrm->avail_out == OUT_CHUNK) return RL_END;
				if (ret < 0) return RL_ERR;

				if (!line_start)
					line_start = (char*)z_out.data();
			}

			//line_end = strstr(line_start, "\n"); // returns 0 if '\n' not found
			line_end = std::find(line_start, (char*)&z_out[0] + OUT_CHUNK - pstrm->avail_out - 1, 10); // '\n'
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_CHUNK - strm.avail_out - 1, [l = std::locale{}](auto ch) { return ch == 10; });
			//line_end = std::find_if(line_start, (char*)&z_out[0] + OUT_CHUNK - strm.avail_out - 1, [l = std::locale{}](auto ch) { return std::isspace(ch, l); });
			if (line_end && line_end[0] == 10)
			{
				std::copy(line_start, line_end, std::back_inserter(line));

				if (line_end < (char*)&z_out[0] + OUT_CHUNK - pstrm->avail_out - 1) // check there is data after line_end
					line_start = line_end + 1; // skip '\n'
				else
					line_start = 0; // no more data in out buffer - flag to inflate more

				line_ready = true;
			}
			else
			{
				line_end = (char*)&z_out[0] + OUT_CHUNK - pstrm->avail_out; // end of data in out buffer
				std::copy(line_start, line_end, std::back_inserter(line));
				line_start = (pstrm->avail_out == 0) ? 0 : line_end;
				line_end = 0;
				line_ready = false;
			}

			//if (line == ">L1S29.274586_481951 HWI-EAS440_0386:1:28:10252:18627#0/1")
			//	std::cout << "HERE\n";
		} // ~for !line_ready
	}
	else
		std::getline(ifs, line);

	return RL_OK;
} // ~Reader::getline

int Reader::inf(std::ifstream & ifs)
{
	int ret;
	std::stringstream ss;

	for (;;)
	{
		if (pstrm->avail_in == 0 && !ifs.eof()) // in buffer empty
		{
			std::fill(z_in.begin(), z_in.end(), 0); // reset buffer to 0
			ifs.read((char*)z_in.data(), CHUNK);
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
			if (pstrm->avail_out < OUT_CHUNK)	pstrm->avail_out = OUT_CHUNK;
			inflateEnd(pstrm);
			return Z_STREAM_END;
		}

		if (pstrm->avail_out == 0) // out buffer is full - reset
		{
			std::fill(z_out.begin(), z_out.end(), 0); // reset buffer to 0
			pstrm->avail_out = OUT_CHUNK;
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

		if (pstrm->avail_out == 0 || (pstrm->avail_out < OUT_CHUNK && pstrm->avail_in == 0)) break;
	} // for(;;)

	return ret;// == Z_STREAM_END ? Z_OK : Z_DATA_ERROR;
} // ~Reader::inf

void Reader::initZstream()
{
	static z_stream strm; // TODO: better way?
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

	z_in.resize(CHUNK);
	z_out.resize(OUT_CHUNK);
	std::fill(z_in.begin(), z_in.end(), 0);
	std::fill(z_out.begin(), z_out.end(), 0);
} // ~Reader::initZstream