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
#include <thread>

#include "readsfile.hpp"
#include "izlib.hpp"
#include "common.hpp"

Readsfile::Readsfile(std::vector<std::string>& readfiles, bool is_gz)
	:
	is_done(false),
	count_all(0),
	is_gzipped(is_gz),
	is_two_files(readfiles.size() > 1),
	is_next_fwd(true),
	readfiles(readfiles),
	izlib_fwd(is_gz),
	izlib_rev(is_gz)
{} // ~Readsfile::Readsfile

//Readsfile::~Readsfile() {}

/* 
 * thread runnable 
 */
void Readsfile::run()
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("Readsfile::run thread ", std::this_thread::get_id(), " started");

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

	// bug 114 - Z_STREAM_ERROR even though init is called at construction time
	izlib_fwd.init();
	izlib_rev.init();

	// loop until EOF - get reads - push on queue
	for (;;)
	{
		if (is_next_fwd) {
			if (nextfwd(fs_fwd).size() > 0) {
				++count_all;
				if (is_two_files)
					is_next_fwd = false;
			}
		}
		else {
			if (nextrev(fs_rev).size() > 0) {
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
} // ~Readsfile::run

bool Readsfile::loadReadByIdx(Read & read)
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

} // ~Readsfile::loadRead

bool Readsfile::loadReadById(Read& read)
{
	return true;
} // ~Readsfile::loadReadById

std::string Readsfile::next(std::ifstream& ifs) {
	std::string line;
	return line;
}

/* 
 * @param array of read files - one or two (if paired) files
 * @return string with format: 'read_id \n read', where read_id = 'filenum_readnum' e.g. '0_1001', read = 'header \n sequence \n quality'
 */
std::string Readsfile::nextfwd(std::ifstream& ifs) {

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
			stat = izlib_fwd.getline(ifs, line);

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
				read << file_num << '_' << state_fwd.read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
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
} // ~Readsfile::nextfwd

/* 
 * TODO: identical to nextfwd.
 * dereferencing Gzip always causes errors. 
 * Cannot store Gzip in an array and cannot pass it by reference (may be can).
 */
std::string Readsfile::nextrev(std::ifstream& ifs) 
{
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
			stat = izlib_rev.getline(ifs, line);

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
				read << file_num << '_' << state_rev.read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
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
} // ~Readsfile::nextrev


/**
   get a next read sequence from a reads file

   @param  fstreams IN     vector of reads streams, one entry per a reads file
   @param  vzlib    IN     vector of Izlib interfaces, one entry per a reads file
   @param  vstate   IN     vector of reading states, one per a reads file
   @param  inext    IN/OUT index of the stream to read. Automatically toggled.
   @param  seq      OUT    read sequence
   @return true if record exists, else false
 */
bool Readsfile::next(std::vector<std::ifstream>& fstreams, 
                     std::vector<Izlib>& vzlib, 
                     std::vector<Readstate>& vstate, 
                     int& inext, 
                     std::string& seq )
{
	std::string line;
	std::stringstream read; // an empty read
	auto stat = vstate[inext].last_stat;

	// read lines from the reads file and extract a single read
	for (auto count = vstate[inext].last_count; !vstate[inext].is_done; ++count) // count lines in a single record/read
	{
		if (vstate[inext].last_header.size() > 0)
		{
			read << inext << '_' << vstate[inext].read_count << '\n';
			read << vstate[inext].last_header << '\n';
			vstate[inext].last_header = "";
		}

		// read a line
		if (!vstate[inext].is_done)
			stat = izlib_rev.getline(fstreams[inext], line);

		// EOF reached - return last read
		if (stat == RL_END)
		{
			vstate[inext].is_done = true;
			auto FR = inext == 0 ? "FWD" : "REF";
			INFO("EOF ", FR, " reached. Total reads: ", ++vstate[inext].read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			auto FR = inext == 0 ? "FWD" : "REF";
			ERR("reading from ", FR, " file. Exiting...");
			exit(1);
		}

		if (line.empty())
		{
			--count;
			continue;
		}

		++vstate[inext].line_count;

		// right-trim whitespace in place (removes '\r' too)
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());

		// the first line in file
		if (vstate[inext].line_count == 1)
		{
			vstate[inext].isFastq = (line[0] == FASTQ_HEADER_START);
			vstate[inext].isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && vstate[inext].isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((vstate[inext].isFasta && line[0] == FASTA_HEADER_START) || (vstate[inext].isFastq && count == 0)) // header line reached
		{
			if (vstate[inext].line_count == 1)
			{
				read << inext << '_' << vstate[inext].read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
				read << line << '\n'; // the very first header
				count = 0;
			}
			else
			{
				// read is ready - return
				vstate[inext].last_header = line;
				vstate[inext].last_count = 1;
				vstate[inext].last_stat = stat;
				break;
			}
		} // ~if header line
		else
		{ // add sequence -->
			if (vstate[inext].isFastq)
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

	++vstate[inext].read_count;

	if (fstreams.size() > 1)
		inext = inext == 0 ? 1 : 0; // toggle the stream index

	return read.str().size() > 0;
} // ~Readsfile::next

/**
 * test if there is a next read in the reads file
 */
bool Readsfile::hasnext(std::ifstream& ifs)
{
	bool is_next = false;
	return is_next;
} // ~Readsfile::hasnext

void Readsfile::reset()
{
	is_done = false;
	count_all = 0;
	is_next_fwd = true;
}

/*
 input options:
   2 paired file
   1 paired file
   1 non-paired file
*/
bool Readsfile::split(const unsigned num_parts, const std::string& outdir)
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("start splitting");

	std::vector<std::ifstream> fstreams(readfiles.size());
	for (int i = 0; i < readfiles.size(); ++i) {
		fstreams[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
		if (!fstreams[i].is_open()) {
			std::cerr << STAMP << "Failed to open file " << readfiles[i] << std::endl;
			exit(1);
		}
	}

	std::vector<Izlib> vzlib(2, Izlib(true));
	std::vector<Readstate> vstate(2);
	int inext = 0; // fwd
	std::string readstr;

	// loop until EOF - get reads - push on queue
	for (;;)
	{
		if (next(fstreams, vzlib, vstate, inext, readstr)) {
			++count_all;
			if (is_two_files)
				is_next_fwd = false;
		}

		is_done = true;
		for (auto state: vstate ) {
			is_done = is_done && state.is_done;
		}
		if (is_done) {
			break;
		}
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("Done splitting. Reads count: ", count_all, " Runtime sec: ", elapsed.count(), "\n");
	return true;
} // ~Readsfile::split