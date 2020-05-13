/**
 * FILE: reader.cpp
 * Created: Nov 26, 2017 Sun
 * @copyright 2016-20 Clarity Genomics BVBA
 * 
 * Processes Reads file, creates Read objects, and pushes them to a queue for further pick−up by Processor
 */

#include <vector>
#include <iostream>
#include <sstream> // std::stringstream
#include <ios> // std::ios_base
#include <algorithm> // find, find_if
#include <chrono> // std::chrono
#include <iomanip> // std::precision
#include <locale> // std::isspace

#include "reader.hpp"
#include "gzip.hpp"

Reader::Reader(ReadsQueue& readQueue, std::vector<std::string>& readfiles, bool is_gz)
	:
	is_done(false),
	count_all(0),
	is_gzipped(is_gz),
	next_idx(0),
	readfiles(readfiles),
	readQueue(readQueue)
{
	for (auto i = 0; i < readfiles.size(); ++i) {
		states.reserve(readfiles.size());
		states.emplace_back(Readstate(is_gz));
	}
} // ~Reader::Reader

//Reader::~Reader() {}

/* 
 * thread runnable 
 */
void Reader::run()
{
	{
		std::stringstream ss;
		ss << STAMP << "Reader::run " << "thread " << std::this_thread::get_id() << " started" << std::endl;
		std::cout << ss.str();
	}

	// open file streams
	std::vector<std::ifstream> fsl(readfiles.size());
	for (auto i = 0; i < readfiles.size(); ++i) {
		fsl[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
		if (!fsl[i].is_open()) {
			std::cerr << STAMP << "Failed to open file " << readfiles[i] << std::endl;
			exit(1);
		}
	}

	// loop until EOF - get reads - push on queue
	for (bool is_ok = false; !is_done;)
	{
		if (is_done) {
			readQueue.is_done_push = true;
			{
				std::stringstream ss;
				ss << STAMP << "Reader::run " << " Done Reading from all streams" << std::endl;
				std::cout << ss.str();
			}
			break;
		}
		is_ok = readQueue.push(nextread(fsl));
		if (is_ok) ++count_all;
	}

	{
		std::stringstream ss;
		ss << "Reader::run " << " thread " << std::this_thread::get_id() << " done. Reads count: " << count_all << std::endl;
		std::cout << ss.str();
	}
} // ~Reader::run

bool Reader::loadReadByIdx(Read & read)
{
	std::stringstream ss;
	bool isok = false;
#if 0

	pstate->fs.open(Reader::*pstate->readsfile, std::ios_base::in | std::ios_base::binary);
	if (!rstate->fs.is_open()) {
		std::cerr << STAMP << "failed to open " << rstate->readsfile << std::endl;
		exit(EXIT_FAILURE);
	}

	std::string line;
	std::size_t read_num = 0;
	//bool isFastq = true;
	//Gzip gzip(opts.is_gz);

	auto t = std::chrono::high_resolution_clock::now();

	// read lines from the reads file
	for (int count = 0, stat = 0; ; ) // count lines in a single read
	{
		stat = rstate->gzip.getline(rstate->fs, line);
		if (stat == RL_END) break;

		if (stat == RL_ERR)
		{
			std::cerr << STAMP << "ERROR reading from Reads file. Exiting..." << std::endl;
			exit(1);
		}

		if (line.empty()) continue;

		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());

		if ( line[0] == FASTA_HEADER_START || line[0] == FASTQ_HEADER_START )
		{
			if (!read.isEmpty) {
				isok = true;
				break; // read is ready
			}

			// add header -->
			if (read_num == read.read_num)
			{
				rstate->isFastq = (line[0] == FASTQ_HEADER_START);
				read.format = rstate->isFastq ? Format::FASTQ : Format::FASTA;
				read.header = line;
				read.isEmpty = false;
			}
			else {
				++read_num;
				count = 0; // for fastq
			}
		} // ~if header line
		else if ( !read.isEmpty )
		{
			// add sequence -->
			if (rstate->isFastq )
			{
				++count;
				if ( line[0] == '+' ) continue;
				if ( count == 3 )
				{
					read.quality = line; // last line in Fastq read
					continue;
				}
			}
			read.sequence += line;
		}
	} // ~for getline

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;

	//ss << id << " thread: " << std::this_thread::get_id() << " done. Elapsed time: "
	//	<< std::setprecision(2) << std::fixed << elapsed.count() << " sec Reads added: " << read_id << std::endl;
	//std::cout << ss.str(); ss.str("");

	rstate->fs.close();
#endif
	return isok;

} // ~Reader::loadRead

bool Reader::loadReadById(Read & read)
{
	return true;
} // ~Reader::loadReadById

/* 
 * @param array of read files - one or two (if paired) files
 * @return string with format: 'read_id \n read', where read_id = 'filenum_readnum' e.g. '0_1001', read = 'header \n sequence \n quality'
 */
std::string Reader::nextread(std::vector<std::ifstream>& fsl) {

	std::string line;
	std::stringstream read; // an empty read
	auto stat = states[next_idx].last_stat;

	// read lines from the reads file and extract a single read
	for (auto count = states[next_idx].last_count; !states[next_idx].is_done; ++count) // count lines in a single record/read
	{
		if (states[next_idx].last_header.size() > 0)
		{
			read << next_idx << '_' << states[next_idx].read_count << '\n';
			read << states[next_idx].last_header << '\n';
			states[next_idx].last_header = "";
		}

		// read a line
		if (!states[next_idx].is_done)
			stat = states[next_idx].gzip.getline(fsl[next_idx], line);


		// EOF reached - return last read
		if (stat == RL_END)
		{
			states[next_idx].is_done = true;
			is_done = readfiles.size() == 2 ? states[0].is_done & states[1].is_done : states[0].is_done;
			{
				std::stringstream ss;
				ss << STAMP << "EOF reached. File index: " << next_idx << " Total reads: " << ++states[next_idx].read_count << std::endl;
				std::cout << ss.str();
			}
			break;
		}

		if (stat == RL_ERR)
		{
			std::cerr << STAMP << "ERROR reading from file: [" << readfiles[next_idx] << "]. Exiting..." << std::endl;
			exit(1);
		}

		if (line.empty())
		{
			--count;
			continue;
		}

		++states[next_idx].line_count;

		// right-trim whitespace in place (removes '\r' too)
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());

		// the first line in file
		if (states[next_idx].line_count == 1)
		{
			states[next_idx].isFastq = (line[0] == FASTQ_HEADER_START);
			states[next_idx].isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && states[next_idx].isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((states[next_idx].isFasta && line[0] == FASTA_HEADER_START) || (states[next_idx].isFastq && count == 0)) // header line reached
		{
			if (states[next_idx].line_count == 1)
			{
				read << next_idx << '_' << states[next_idx].read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
				read << line << '\n'; // the very first header
				count = 0;
			}
			else 
			{
				// read is ready - return
				states[next_idx].last_header = line;
				states[next_idx].last_count = 1;
				states[next_idx].last_stat = stat;
				break;
			}
		} // ~if header line
		else
		{ // add sequence -->
			if (states[next_idx].isFastq)
			{
				if (count == 2) // line[0] == '+' validation is already done by readstats::calculate
					continue;
				if (count == 3)
				{
					read << line;
					continue;
				}
				read << line << '\n'; // FQ sequence
			}
			else {
				read << line; // FASTA sequence possibly multiline
			}
		}
	} // ~for getline

	++states[next_idx].read_count;

	// toggle next file
	if (readfiles.size() == 2) {
		next_idx = next_idx == 0 ? 1 : 0;
	}

	return read.str();
} // ~Reader::nextread


/**
 * get a next read sequence from the reads file
 * @return true if record exists, else false
 */
bool Reader::nextread(std::string readsfile, std::string &seq)
{
	bool has_seq = false;
	seq = ""; // ensure empty
	std::string line;
#if 0
	// read lines from the reads file and create Read object
	for (int count = last_count, stat = last_stat; !is_done; ++count) // count lines in a single record/read
	{
		stat = gzip.getline(ifs, line);

		if (stat == RL_ERR)
		{
			std::cerr << STAMP << "ERROR reading from file: [" << readsfile << "]. Exiting..." << std::endl;
			exit(1);
		}

		if (stat == RL_END)
		{
			// the last Read processed
			if (!seq.empty())
			{
				++read_count;
				is_done = true;
				reset();
			}
			break;
		}

		// skip empty lines
		if (line.empty())
		{
			--count;
			continue;
		}

		++line_count; // total lines read so far

		// left trim space and '>' or '@'
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());

		if (line_count == 1)
		{
			isFastq = (line[0] == FASTQ_HEADER_START);
			isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((isFasta && line[0] == FASTA_HEADER_START) || (isFastq && count == 0)) // header reached
		{
			if (!seq.empty())
			{
				++read_count;
				last_count = 1;
				last_stat = stat;
				break;
			}

			count = 0; // FASTA record start
		} // ~if header line
		else
		{
			if (isFastq && (count == 2 || count == 3))
			{
				continue; // skip + and quality line in Fastq
			}

			seq += line; // FASTA multi-line sequence or FASTQ sequence
		}
	} // ~for getline

	has_seq = seq.size() > 0 ? true : false;
#endif
	return has_seq;
} // ~Reader::nextread

/**
 * test if there is a next read in the reads file
 */
bool Reader::hasnext(std::ifstream& ifs)
{
	bool is_next = false;
	return is_next;
} // ~Reader::hasnext

void Reader::reset()
{
	is_done = false;
	next_idx = 0;
}