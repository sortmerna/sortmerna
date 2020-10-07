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
#include <filesystem>

#include "readfeed.hpp"
#include "izlib.hpp"
#include "common.hpp"

/*
 @param type       feed type
 @param readfiles  vector with reads file paths
 @param is_gz      flags the readsfiles format gzip | non-gzip (flat)
*/
Readfeed::Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, bool is_gz)
	:
	type(type),
	is_done(false),
	is_ready(false),
	count_all(0),
	is_gzipped(is_gz),
	is_two_files(readfiles.size() > 1),
	is_next_fwd(true),
	readfiles(readfiles)
{
	// verify the read feed is ready
	if (type == FEED_TYPE::SPLIT_READS)
		is_ready = is_split_ready();
} // ~Readfeed::Readfeed 1

Readfeed::Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, const unsigned num_parts, const unsigned num_reads, const std::string& outdir)
	:
	type(type),
	is_done(false),
	is_ready(false),
	count_all(0),
	is_gzipped(false),
	is_two_files(readfiles.size() > 1),
	is_next_fwd(true),
	readfiles(readfiles)
{
	// define files format (FASTA, FASTQ) and compression (gz, non-gz)

} //~Readfeed::Readfeed 2

//Readfeed::~Readfeed() {}

/* 
  verify the split was already performed and the feed is ready
  Split readfeed descriptor:
    timestamp: xxx
	num_input: 2  # number of input files
	num_parts: 3  # number of split parts. split[].size = num_input * num_parts
    input:
	  - file_1: name, sha
	  - file_2: name, sha
	split:
	 - file_1: name, sha
	 - file_2: name, sha
	 ...
	 - file_n: name, sha
*/
bool Readfeed::is_split_ready() {
	is_ready = false;
	return is_ready;
} // ~Readfeed::is_split_ready

/* 
 * can be run in a thread
 */
void Readfeed::run()
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("Readsfile::run thread ", std::this_thread::get_id(), " started");

	// init input file streams
	ifsv.resize(readfiles.size());
	for (int i = 0; i < readfiles.size(); ++i) {
		ifsv[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", readfiles[i]);
			exit(1);
		}
	}

	// bug 114 - Z_STREAM_ERROR even though init is called at construction time
	vzlib_in.resize(readfiles.size(), Izlib(true));
	for (auto i = 0; i < vzlib_in.size(); ++i) {
		vzlib_in[i].init(false);
	}

	std::string readstr;
	int inext = 0;

	// loop until EOF - get reads - push on queue
	for (;;)
	{
		if (next(inext, readstr)) {
			++count_all;
			if (is_two_files) {
				is_next_fwd = false;
				inext = inext == 0 ? 1 : 0; // toggle next file index
			}
		}

		is_done = true;
		for (auto i = 0; i < vstate_in.size(); ++i) {
			is_done = is_done && vstate_in[i].is_done;
		}
		if (is_done) {
			break;
		}
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("Reader::run thread ", std::this_thread::get_id(), " done. Reads count: ", count_all, " Runtime sec: ", elapsed.count());
} // ~Readfeed::run

bool Readfeed::loadReadByIdx(Read & read)
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

} // ~Readfeed::loadRead

bool Readfeed::loadReadById(Read& read)
{
	return true;
} // ~Readfeed::loadReadById

/*
 @param fs_idx  index of the file stream into Readfeed::ifsv
*/
bool Readfeed::next(int fs_idx, std::string& seq) {
	std::string line;
	std::stringstream read; // an empty read
	auto stat = vstate_in[fs_idx].last_stat;
	auto inext = vnext_fwd_in[fs_idx] ? fs_idx : fs_idx + 1; // index of the next file stream i.e. idx fwd | idx rev

	// read lines from the reads file and extract a single read
	for (auto count = vstate_in[inext].last_count; !vstate_in[inext].is_done; ++count) // count lines in a single record/read
	{
		if (vstate_in[inext].last_header.size() > 0)
		{
			read << inext << '_' << vstate_in[inext].read_count << '\n';
			read << vstate_in[inext].last_header << '\n';
			vstate_in[inext].last_header = "";
		}

		// read a line
		if (!vstate_in[inext].is_done)
			stat = vzlib_in[inext].getline(ifsv[inext], line);

		// EOF reached - return last read
		if (stat == RL_END)
		{
			vstate_in[inext].is_done = true;
			auto FR = inext == 0 ? "FWD" : "REF";
			INFO("EOF ", FR, " reached. Total reads: ", ++vstate_in[inext].read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			auto FR = vnext_fwd_in[fs_idx] ? "FWD" : "REF";
			ERR("reading from ", FR, " file. Exiting...");
			exit(1);
		}

		if (line.empty())
		{
			--count;
			continue;
		}

		++vstate_in[inext].line_count;

		// right-trim whitespace in place (removes '\r' too)
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());

		// the first line in file
		if (vstate_in[inext].line_count == 1)
		{
			vstate_in[inext].isFastq = (line[0] == FASTQ_HEADER_START);
			vstate_in[inext].isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && vstate_in[inext].isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((vstate_in[inext].isFasta && line[0] == FASTA_HEADER_START) || (vstate_in[inext].isFastq && count == 0)) // header line reached
		{
			if (vstate_in[inext].line_count == 1)
			{
				read << inext << '_' << vstate_in[inext].read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
				read << line << '\n'; // the very first header
				count = 0;
			}
			else
			{
				// read is ready - return
				vstate_in[inext].last_header = line;
				vstate_in[inext].last_count = 1;
				vstate_in[inext].last_stat = stat;
				break;
			}
		} // ~if header line
		else
		{ // add sequence -->
			if (vstate_in[inext].isFastq)
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

	++vstate_in[inext].read_count;

	vnext_fwd_in[fs_idx] != vnext_fwd_in[fs_idx]; // toggle the stream index

	if (read.str().size() > 0)
		seq = read.str();
	return read.str().size() > 0;
} // ~Readfeed::next


/**
   get a next read sequence from a reads file

   @param  fstreams IN     vector of reads streams, one entry per a reads file
   @param  vzlib    IN     vector of Izlib interfaces, one entry per a reads file
   @param  vstate   IN     vector of reading states, one per a reads file
   @param  inext    IN/OUT index of the stream to read. Automatically toggled.
   @param  seq      OUT    read sequence
   @return true if record exists, else false
 */
bool Readfeed::next(std::vector<std::ifstream>& fstreams,
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
			stat = vzlib[inext].getline(fstreams[inext], line);

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

	if (read.str().size() > 0)
		seq = read.str();
	return read.str().size() > 0;
} // ~Readfeed::next

/**
 * test if there is a next read in the reads file
 */
bool Readfeed::hasnext(std::ifstream& ifs)
{
	bool is_next = false;
	return is_next;
} // ~Readfeed::hasnext

void Readfeed::reset()
{
	is_done = false;
	count_all = 0;
	is_next_fwd = true;
}

/*
 input options:
   2 paired file     -> num files out = 2 x num_parts i.e. each paired file is split into specified number of parts
   1 paired file     -> num files out = 1 x num_parts AND ensure each part has even number of reads
   1 non-paired file -> num files out = 1 x num_parts

  @param num_parts  number of parts to split the file into
  @param num_reads  total number of reads in all read files (Readstats::all_reads_count)
*/
bool Readfeed::split(const unsigned num_parts, const unsigned num_reads, const std::string& outdir)
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("start splitting");
	auto retval = true;

	//std::vector<std::ifstream> ifsv(readfiles.size());
	ifsv.resize(readfiles.size());
	for (int i = 0; i < readfiles.size(); ++i) {
		ifsv[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", readfiles[i]);
			exit(1);
		}
	}

	// prepare zlib interfaces for inflating reads files
	bool is_gz = true;
	//std::vector<Izlib> vzlib_in(2, Izlib(is_gz));
	vzlib_in.resize(readfiles.size(), Izlib(is_gz));
	for (auto i = 0; i < vzlib_in.size(); ++i) {
		vzlib_in[i].init(false);
	}

	//std::vector<Readstate> vstate_in(2);
	vstate_in.resize(readfiles.size());
	int inext = 0; // fwd - index of the input file
	std::string readstr;

	// prepare split files to output
	//std::vector<std::ofstream> ofsv(num_parts * readfiles.size());
	ofsv.resize(num_parts * readfiles.size());
	size_t idx = 0; // stream index
	for (int i = 0; i < readfiles.size(); ++i) {
		auto pdir = std::filesystem::path(readfiles[i]).parent_path();
		auto stem = i == 0 ? "fwd_" : "rev_";
		std::stringstream ss;
		for (int j = 0; j < num_parts; ++j) {
			ss << stem << j << ".fq.gz";
			auto fn = pdir / ss.str();
			ofsv[idx].open(fn, std::ios::app | std::ios::binary);
			if (!ofsv[idx].good()) {
				ERR("Failed to open file ", fn);
				exit(1);
			}
			else {
				INFO("opened file: ", fn);
			}
			ss.str("");
			++idx;
		}
	}

	// prepare zlib interface for writing split files
	//std::vector<Izlib> vzlib_out(num_parts * readfiles.size(), Izlib(true, true, true));
	vzlib_out.resize(num_parts * readfiles.size(), Izlib(true, true, true));
	for (auto i = 0; i < vzlib_out.size(); ++i) {
		vzlib_out[i].init(true);
	}

	//std::vector<Readstate> vstate_out(num_parts * readfiles.size());
	vstate_out.resize(num_parts * readfiles.size());

	// calculate number of reads in each of the output files
	auto nreads = readfiles.size() == 2 ? num_reads / 2 : num_reads; // num reads in a single input file e.g. FWD
	auto minr = nreads / num_parts; // quotient i.e. min number of reads in each output file
	auto surplus = nreads - minr * num_parts; // remainder of reads to be distributed between the output files
	for (auto i = 0; i < num_parts; ++i) {
		auto maxr = i < surplus ? minr + 1 : minr; // distribute the surplus
		for (auto j = 0; j < readfiles.size(); ++j) {
			vstate_out[i + j * num_parts].max_reads = maxr;
		}
	}

	int iout = 0;

	// loop until EOF - get reads - write into split files
	for (;;)
	{
		if (next(ifsv, vzlib_in, vstate_in, inext, readstr)) {
			++count_all;
			if (is_two_files) {
				is_next_fwd = false;
				++vstate_out[iout].read_count;
				auto ret = vzlib_out[iout].defstr(readstr, ofsv[iout], vstate_out[iout].read_count == vstate_out[iout].max_reads); // Z_STREAM_END | Z_OK - ok
				if (ret < Z_OK || ret > Z_STREAM_END) {
					ERR("Failed deflating readstring: ", readstr, " Output file idx: ", iout, " zlib status: ", ret);
					retval = false;
					break;
				}
				// set next value of the out file index
				if (iout == 2 * num_parts - 1) iout = 0;
				else {
					iout = inext == 0 ? iout - num_parts + 1 : iout + num_parts;
				}
			}
		}

		is_done = true;
		for (auto i = 0; i < vstate_in.size(); ++i) {
			is_done = is_done && vstate_in[i].is_done;
		}
		if (is_done) {
			break;
		}
	} // ~for

	// close in file streams
	for (int i = 0; i < ifsv.size(); ++i) {
		if (ifsv[i].is_open()) {
			ifsv[i].close();
		}
	}

	// close out file streams
	for (auto i = 0; i < ofsv.size(); ++i) {
		if (ofsv[i].is_open())
			ofsv[i].close();
	}

	// clean up zlib deflate streams
	for (auto i = 0; i < vzlib_out.size(); ++i) {
		vzlib_out[i].clean();
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("Done splitting. Reads count: ", count_all, " Runtime sec: ", elapsed.count(), "\n");
	return retval;
} // ~Readfeed::split

/*
  20200924 Thu
  verify the split is already performed and ready to use

  Upon the split the following descriptor file is generated in the split files directory:
    stamp:
	  time: <timestamp>
      source files count: 2
      split files count:  6
	  total reads:        10M
	  source read files: []
	  split files: []

  @param  num_part   number of splits on each source reads file
  @param  num_reads  total count of reads in all source files
  @param  dbdir      directory where the split reads are located
  @param  readfiles  source reads files paths
*/
bool Readfeed::is_split_done(
	const unsigned num_parts, 
	const unsigned num_reads, 
	const std::string dbdir, 
	const std::vector<std::string>& readfiles)
{
	INFO("TODO: implement");
	return true;
} // ~Readfeed::is_split_ready