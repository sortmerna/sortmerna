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
	is_two_files(readfiles.size() > 1),
	is_next_fwd(true),
	readfiles(readfiles),
	readQueue(readQueue),
	gzip_fwd(is_gz),
	gzip_rev(is_gz)
{} // ~Reader::Reader

//Reader::~Reader() {}

/* 
 * thread runnable 
 */
void Reader::run()
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("Reader::run thread ", std::this_thread::get_id(), " started");

	std::ifstream fs_fwd;
	std::ifstream fs_rev;
	fs_fwd.open(readfiles[0], std::ios_base::in | std::ios_base::binary);
	if (!fs_fwd.is_open()) {
		std::cerr << STAMP << "Failed to open file " << readfiles[0] << std::endl;
		exit(1);
	}
	if (is_two_files) {
		fs_rev.open(readfiles[1], std::ios_base::in | std::ios_base::binary);
		if (!fs_rev.is_open()) {
			std::cerr << STAMP << "Failed to open file " << readfiles[1] << std::endl;
			exit(1);
		}
	}

	// without this getting Z_STREAM_ERROR even though init is called at construction time
	gzip_fwd.init();
	gzip_rev.init();

	// loop until EOF - get reads - push on queue
	for (;;)
	{
		if (is_next_fwd) {
			if (readQueue.push(nextfwd(fs_fwd))) {
				++count_all;
				if (is_two_files)
					is_next_fwd = false;
			}
		}
		else {
			if (readQueue.push(nextrev(fs_rev))) {
				++count_all;
				is_next_fwd = true;
			}
		}
		
		is_done = is_two_files ? state_fwd.is_done && state_rev.is_done : state_fwd.is_done;
		if (is_done) {
			INFO("Reader::run thread ", std::this_thread::get_id(), " Done Reading from all streams");
			break;
		}
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("Reader::run thread ", std::this_thread::get_id(), " done. Reads count: ", count_all, " Runtime sec: ", elapsed.count());
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

bool Reader::loadReadById(Read& read)
{
	return true;
} // ~Reader::loadReadById

std::string Reader::nextread(std::ifstream& ifs) {
	std::string line;
	return line;
}

/* 
 * @param array of read files - one or two (if paired) files
 * @return string with format: 'read_id \n read', where read_id = 'filenum_readnum' e.g. '0_1001', read = 'header \n sequence \n quality'
 */
std::string Reader::nextfwd(std::ifstream& ifs) {

	std::string line;
	std::stringstream read; // an empty read
	auto stat = state_fwd.last_stat;
	auto file_num = 0; // FWD file

	// read lines from the reads file and extract a single read
	for (auto count = state_fwd.last_count; !state_fwd.is_done; ++count) // count lines in a single record/read
	{
		if (state_fwd.last_header.size() > 0)
		{
			read << file_num << '_' << state_fwd.read_count << '\n';
			read << state_fwd.last_header << '\n';
			state_fwd.last_header = "";
		}

		// read a line
		if (!state_fwd.is_done)
			stat = gzip_fwd.getline(ifs, line);

		// EOF reached - return last read
		if (stat == RL_END)
		{
			state_fwd.is_done = true;
			INFO("EOF FWD reached. Total reads: ", ++state_fwd.read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			ERR("reading from FWD file. Exiting...");
			exit(1);
		}

		if (line.empty())
		{
			--count;
			continue;
		}

		++state_fwd.line_count;

		// right-trim whitespace in place (removes '\r' too)
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());

		// the first line in file
		if (state_fwd.line_count == 1)
		{
			state_fwd.isFastq = (line[0] == FASTQ_HEADER_START);
			state_fwd.isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && state_fwd.isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((state_fwd.isFasta && line[0] == FASTA_HEADER_START) || (state_fwd.isFastq && count == 0)) // header line reached
		{
			if (state_fwd.line_count == 1)
			{
				read << 0 << '_' << state_fwd.read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
				read << line << '\n'; // the very first header
				count = 0;
			}
			else 
			{
				// read is ready - return
				state_fwd.last_header = line;
				state_fwd.last_count = 1;
				state_fwd.last_stat = stat;
				break;
			}
		} // ~if header line
		else
		{ // add sequence -->
			if (state_fwd.isFastq)
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

	++state_fwd.read_count;

	return read.str();
} // ~Reader::nextfwd

/* 
 * TODO: identical to nextfwd.
 * dereferencing Gzip always causes errors. 
 * Cannot store Gzip in an array and cannot pass it by reference (may be can).
 */
std::string Reader::nextrev(std::ifstream& ifs) {

	std::string line;
	std::stringstream read; // an empty read
	auto stat = state_rev.last_stat;
	auto file_num = 1; // REV file

	// read lines from the reads file and extract a single read
	for (auto count = state_rev.last_count; !state_rev.is_done; ++count) // count lines in a single record/read
	{
		if (state_rev.last_header.size() > 0)
		{
			read << file_num << '_' << state_rev.read_count << '\n';
			read << state_rev.last_header << '\n';
			state_rev.last_header = "";
		}

		// read a line
		if (!state_rev.is_done)
			stat = gzip_rev.getline(ifs, line);

		// EOF reached - return last read
		if (stat == RL_END)
		{
			state_rev.is_done = true;
			INFO("EOF REV reached. Total reads: " , ++state_rev.read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			ERR("reading from REV file. Exiting...");
			exit(1);
		}

		if (line.empty())
		{
			--count;
			continue;
		}

		++state_rev.line_count;

		// right-trim whitespace in place (removes '\r' too)
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());

		// the first line in file
		if (state_rev.line_count == 1)
		{
			state_rev.isFastq = (line[0] == FASTQ_HEADER_START);
			state_rev.isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && state_rev.isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((state_rev.isFasta && line[0] == FASTA_HEADER_START) || (state_rev.isFastq && count == 0)) // header line reached
		{
			if (state_rev.line_count == 1)
			{
				read << 0 << '_' << state_rev.read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
				read << line << '\n'; // the very first header
				count = 0;
			}
			else
			{
				// read is ready - return
				state_rev.last_header = line;
				state_rev.last_count = 1;
				state_rev.last_stat = stat;
				break;
			}
		} // ~if header line
		else
		{ // add sequence -->
			if (state_rev.isFastq)
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

	++state_rev.read_count;

	return read.str();
} // ~Reader::nextrev


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
	count_all = 0;
	is_next_fwd = true;
}