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

 Split reads logic:
   - check the split was already done:
     - check readb directory for the descriptor and file
*/
Readfeed::Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, const std::string& basedir)
	:
	type(type),
	biof(BIO_FORMAT::FASTQ),
	zipf(ZIP_FORMAT::GZIP),
	is_done(false),
	is_ready(false),
	is_format_defined(false),
	num_orig_files(readfiles.size()),
	num_splits(0),
	num_reads_tot(0),
	length_all(0),
	basedir(basedir),
	is_two_files(readfiles.size() > 1),
	is_next_fwd(true),
	readfiles(readfiles)
{
	if (type == FEED_TYPE::SPLIT_READS)
		is_ready = is_split_ready();
	else
		is_ready = true;

	if (!is_ready) {
		define_format();
		count_reads();
	}
} // ~Readfeed::Readfeed 1

Readfeed::Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, const unsigned num_parts, const std::string& basedir)
	:
	type(type),
	biof(BIO_FORMAT::FASTQ),
	zipf(ZIP_FORMAT::GZIP),
	is_done(false),
	is_ready(false),
	is_format_defined(false),
	num_orig_files(readfiles.size()),
	num_splits(num_parts),
	num_reads_tot(0),
	length_all(0),
	basedir(basedir),
	is_two_files(readfiles.size() > 1),
	is_next_fwd(true),
	readfiles(readfiles)
{
	if (type == FEED_TYPE::SPLIT_READS)
		is_ready = is_split_ready();
	else
		is_ready = true;

	if (!is_ready) {
		define_format();
		count_reads();
	}
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
	unsigned seqlen = 0;

	// loop until EOF - get reads - push on queue
	for (;;)
	{
		if (next(inext, seqlen, readstr)) {
			++num_reads_tot;
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
	INFO("Reader::run thread ", std::this_thread::get_id(), " done. Reads count: ", num_reads_tot, " Runtime sec: ", elapsed.count());
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
 uses zlib if files are gzipped or reads flat files

 @param inext   IN   index of the file into the Readfeed::readfiles vector
 @param seqlen  OUT  length of the read sequence
 @param read    OUT  read
*/
bool Readfeed::next(int inext, unsigned& seqlen, std::string& read) {
	std::string line;
	std::stringstream readss; // an empty read
	auto stat = vstate_in[inext].last_stat;
	seqlen = 0;

	// read lines from the reads file and extract a single read
	for (auto count = vstate_in[inext].last_count; !vstate_in[inext].is_done; ++count) // count lines in a single record/read
	{
		if (vstate_in[inext].last_header.size() > 0)
		{
			readss << inext << '_' << vstate_in[inext].read_count << '\n';
			readss << vstate_in[inext].last_header << '\n';
			vstate_in[inext].last_header = "";
		}

		// read a line
		if (!vstate_in[inext].is_done) {
			if (vstate_in[inext].isZip) {
				stat = vzlib_in[inext].getline(ifsv[inext], line);
			}
			else {
				if (ifsv[inext].eof()) 
					stat = RL_END;
				else 
					std::getline(ifsv[inext], line);
			}
		}

		// EOF reached - return last read
		if (stat == RL_END)
		{
			vstate_in[inext].is_done = true;
			auto FR = inext == 0 ? FWD : REV;
			INFO("EOF ", FR, " reached. Total reads: ", ++vstate_in[inext].read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			auto FR = inext == 0 ? FWD : REV;
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
				readss << inext << '_' << vstate_in[inext].read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
				readss << line << '\n'; // the very first header
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
					readss << line;
					continue;
				}
				readss << line << '\n'; // FQ sequence
				seqlen = line.length();
			}
			else {
				readss << line; // FASTA sequence possibly multiline
				seqlen += line.length();
			}
		}
	} // ~for getline

	++vstate_in[inext].read_count;

	//vnext_fwd_in[inext] != vnext_fwd_in[inext]; // toggle the stream index

	if (readss.str().size() > 0)
		read = readss.str();
	return readss.str().size() > 0;
} // ~Readfeed::next


/**
   get a next read string from a reads file
   20201012: TODO: exactly the same as next(int, str, str) but doesn't count sequence length.
                   Counting sequence length could be a very small overhead -> no need for this oveload?

   @param  inext    IN     index of the stream to read.
   @param  readstr  OUT    read sequence
   @return true if record exists, else false
 */
bool Readfeed::next(int inext, std::string& readstr )
{
	std::string line;
	std::stringstream readss; // an empty read
	auto stat = vstate_in[inext].last_stat;

	// read lines from the reads file and extract a single read
	for (auto count = vstate_in[inext].last_count; !vstate_in[inext].is_done; ++count) // count lines in a single record/read
	{
		if (vstate_in[inext].last_header.size() > 0)
		{
			readss << inext << '_' << vstate_in[inext].read_count << '\n';
			readss << vstate_in[inext].last_header << '\n';
			vstate_in[inext].last_header = "";
		}

		// read a line
		if (!vstate_in[inext].is_done) {
			if (vstate_in[inext].isZip) {
				stat = vzlib_in[inext].getline(ifsv[inext], line);
			}
			else {
				if (ifsv[inext].eof())
					stat = RL_END;
				else
					std::getline(ifsv[inext], line);
			}
		}

		// EOF reached - return last read
		if (stat == RL_END)
		{
			vstate_in[inext].is_done = true;
			auto FR = inext == 0 ? FWD : REV;
			INFO("EOF ", FR, " reached. Total reads: ", ++vstate_in[inext].read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			auto FR = inext == 0 ? FWD : REV;
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
				readss << inext << '_' << vstate_in[inext].read_count << '\n'; // add read id 'filenum_readnum' starting with '0_0'
				readss << line << '\n'; // the very first header
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
					readss << line;
					continue;
				}
				readss << line << '\n'; // FQ sequence
			}
			else {
				readss << line; // FASTA sequence possibly multiline
			}
		}
	} // ~for getline

	++vstate_in[inext].read_count;

	//vnext_fwd_in[inext] != vnext_fwd_in[inext]; // toggle the stream index

	auto is_read_ok = readss.str().size() > 0;
	if (is_read_ok)
		readstr = readss.str();
	return is_read_ok;
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
	num_reads_tot = 0;
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
bool Readfeed::split(const unsigned num_parts)
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("start splitting");
	auto retval = true;
	num_splits = num_parts;
	auto num_split_files = num_splits * num_orig_files;
	auto tot_files = num_split_files + num_orig_files; // splits + original read files
	ifsv.reserve(tot_files);

	ifsv.resize(num_orig_files);
	for (int i = 0; i < num_orig_files; ++i) {
		ifsv[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", readfiles[i]);
			exit(1);
		}
	}

	// prepare zlib interfaces to inflate reads files
	bool is_gz = true;
	vzlib_in.resize(num_orig_files, Izlib(is_gz));
	for (auto i = 0; i < vzlib_in.size(); ++i) {
		vzlib_in[i].init(false);
	}

	vstate_in.resize(num_orig_files);
	int inext = 0; // fwd - index of the input file
	std::string readstr;

	// prepare split files to output
	// fwd_1.fq.gz, rev_1.fq.gz; fwd_2.fq.gz, rev_2.fq.gz; ...; fwd_n.fq.gz, ref_n.fq.gz
	ofsv.resize(num_split_files);
	size_t idx = 0; // stream index
	for (int i = 0; i < num_orig_files; ++i) {
		std::string stem = i == 0 ? "fwd_" : "rev_";
		std::string sfx_1 = vstate_in[i].isFasta ? ".fa" : ".fq";
		auto sfx = vstate_in[i].isZip ? sfx_1 + ".gz" : sfx_1;
		std::stringstream ss;
		for (int j = 0; j < num_parts; ++j) {
			ss << stem << j << sfx; // split file basename
			auto fn = std::filesystem::path(basedir) / ss.str(); // split file name
			ofsv[idx].open(fn, std::ios::app | std::ios::binary);
			if (!ofsv[idx].good()) {
				ERR("Failed to open file ", fn);
				exit(1);
			}
			else {
				INFO("opened file: ", fn.generic_string());
				readfiles.emplace_back(fn.generic_string()); // add split file for further processing
			}
			ss.str("");
			++idx;
		}
	}

	// prepare zlib interface for writing split files
	vzlib_out.resize(num_split_files, Izlib(true, true));
	for (auto i = 0; i < vzlib_out.size(); ++i) {
		vzlib_out[i].init(true);
	}

	vstate_out.resize(num_split_files);

	// calculate number of reads in each of the output files
	auto nreads = num_orig_files == 2 ? num_reads_tot / 2 : num_reads_tot; // num reads in a single input file e.g. FWD
	auto minr = nreads / num_parts; // quotient i.e. min number of reads in each output file
	auto surplus = nreads - minr * num_parts; // remainder of reads to be distributed between the output files
	for (auto i = 0; i < num_parts; ++i) {
		auto maxr = i < surplus ? minr + 1 : minr; // distribute the surplus
		for (auto j = 0; j < num_orig_files; ++j) {
			vstate_out[i + j * num_parts].max_reads = maxr;
		}
	}

	int iout = 0;
	unsigned seqlen = 0;

	// loop until EOF - get reads - write into split files
	for (;;)
	{
		if (next(inext, seqlen, readstr)) {
			++num_reads_tot;
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
		vzlib_out[i].reset_deflate();
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("Done splitting. Reads count: ", num_reads_tot, " Runtime sec: ", elapsed.count(), "\n");
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

/*
  - define input files format (FASTA, FASTQ) and compression (gz, non-gz)
  - Count reads
  - Count length of all reads

  logic:
    gz:
	  1F 8B 08 08 61 78 C5 5E 00 03 53 52 52 31 36 33 35 38 36 34 5F 31 5F 35 4B 2E 66 61 73 74 71 00
	   0 byte at 8th position_|  |  |_file name starts in ASCII                   file name end_|  |_ 0 byte
	                    ETX byte_|
	flat:
	  first 100 bytes are ascii (<=127 x7F), first char is '@' (x40), and can infer fasta or fastq
*/
bool Readfeed::define_format()
{
	bool is_ascii = true;
	is_format_defined = true;

	if (ifsv.size() < num_orig_files) {
		ifsv.resize(num_orig_files);
	}
	if (vstate_in.size() < num_orig_files) {
		vstate_in.resize(num_orig_files);
	}
	std::string str(100, '\0');
	unsigned seqlen = 0;
	for (int i = 0; i < num_orig_files; ++i) {
		if (!ifsv[i].is_open()) {
			ifsv[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
		}
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", readfiles[i]);
			exit(1);
		}
		ifsv[i].read(&str[0], 100); // get 100 bytes from the stream
		for (auto i = 0; i < str.size(); ++i) {
			// 20201008 TODO: this is quite adhoc - need a better validation like evalutating the gz, zlib header
			if (str[i] > 127 || str[i] < 0) {
				is_ascii = false; // if any of the char is not ascii -> file is not ascii
				break;
			}
		}
		vstate_in[i].isZip = !is_ascii;
		if (!is_ascii) {
			// init izlib for inflation
			if (vzlib_in.size() < i+1) {
				vzlib_in.emplace_back(Izlib());
			}
			vzlib_in[i].init(false);
		}
		
		ifsv[i].seekg(0); // rewind stream to the start
		// get a read
		if (next(i, seqlen, str)) {
			std::string fmt = vstate_in[i].isZip ? "gzipped" : "flat ASCII";
			if (vstate_in[i].isFasta) {
				biof = BIO_FORMAT::FASTA;
				INFO("file: ", readfiles[i], " is FASTA ", fmt);
			}
			else if (vstate_in[i].isFastq) {
				biof = BIO_FORMAT::FASTQ;
				INFO("file: ", readfiles[i], " is FASTQ ", fmt);
			}
			else {
				is_format_defined = false;
				break;
			}
		}
		else {
			is_format_defined = false;
			break;
		}

		if (!is_format_defined) {
			ERR("Cannot define format for file: ", readfiles[i]);
			exit(1);
		}

		// reset Izlib if neessary
		if (vstate_in[i].isZip) {
			auto stat = vzlib_in[i].reset_inflate();
		}
	}

	return is_format_defined;
} // ~Readfeed::define_format

/*
*/
void Readfeed::count_reads()
{
	auto start = std::chrono::high_resolution_clock::now();
	INFO("started count  ...   ");

	if (!is_format_defined) {
		define_format(); // exits if cannot define
	}

	if (ifsv.size() < num_orig_files) {
		ifsv.resize(num_orig_files);
	}
	if (vstate_in.size() < num_orig_files) {
		vstate_in.resize(num_orig_files);
	}
	else {
		for (auto i = 0; i < num_orig_files; ++i) {
			vstate_in[i].reset(); // reset read states
		}
	}

	// prepare file streams
	for (int i = 0; i < num_orig_files; ++i) {
		if (!ifsv[i].is_open()) {
			ifsv[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
		}
		if (ifsv[i].is_open()) {
			ifsv[i].seekg(0); // rewind to the start
		}
		else {
			ERR("Failed to open file ", readfiles[i]);
			exit(1);
		}
	}

	// prepare zlib
	for (int i = 0; i < num_orig_files; ++i) {
		if (vstate_in[i].isZip) {
			vzlib_in[i].init();
		}
	}

	std::string readstr;
	unsigned seqlen = 0;

	// loop until EOF - count reads and sequence lengths
	for (int inext = 0; next(inext, seqlen, readstr);)
	{
		++num_reads_tot;
		length_all += seqlen;
		inext = is_two_files ? inext ^ 1 : inext; // toggle the index of the input file
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("done count. Elapsed time: ", elapsed.count(), " sec. Total reads: ", num_reads_tot);

} // ~Readfeed::count_reads