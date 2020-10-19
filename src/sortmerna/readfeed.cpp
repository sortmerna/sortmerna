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
Readfeed::Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, std::filesystem::path& basedir)
	:
	type(type),
	biof(BIO_FORMAT::FASTQ),
	zipf(ZIP_FORMAT::GZIP),
	is_done(false),
	is_ready(false),
	is_format_defined(false),
	is_two_files(readfiles.size() > 1),
	num_orig_files(readfiles.size()),
	num_splits(0),
	num_reads_tot(0),
	length_all(0),
	min_read_len(0),
	max_read_len(0),
	basedir(basedir),
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

Readfeed::Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, const unsigned num_parts, std::filesystem::path& basedir)
	:
	type(type),
	biof(BIO_FORMAT::FASTQ),
	zipf(ZIP_FORMAT::GZIP),
	is_done(false),
	is_ready(false),
	is_format_defined(false),
	is_two_files(readfiles.size() > 1),
	num_orig_files(readfiles.size()),
	num_splits(num_parts),
	num_reads_tot(0),
	length_all(0),
	min_read_len(0),
	max_read_len(0),
	basedir(basedir),
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
			inext = is_two_files ? inext ^ 1 : inext; // toggle next file index
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
 @param read    OUT  read string
 @param is_orig IN   flags to return the original read string. If false, then the read string has format: 'read_id \n header \n sequence [\n quality]'
*/
bool Readfeed::next(int inext, unsigned& seqlen, std::string& read, bool is_orig) {
	std::string line;
	std::stringstream readss; // an empty read
	auto stat = vstate_in[inext].last_stat;
	seqlen = 0;

	// read lines from the reads file and extract a single read
	for (auto count = vstate_in[inext].last_count; !vstate_in[inext].is_done; ++count) // count lines in a single record/read
	{
		if (vstate_in[inext].last_header.size() > 0)
		{
			if (!is_orig)
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
				if (!is_orig)
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
				if (count == 2) {
					if (is_orig)
						readss << line << '\n'; // line[0] == '+'
					continue;
				}
				if (count == 3)
				{
					readss << line; // FQ quality
					if (is_orig)
						readss << '\n';
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

	if (readss.str().size() > 0)
		read = readss.str();
	return readss.str().size() > 0;
} // ~Readfeed::next


/**
   get a next read string from a reads file
   20201012: TODO: exactly the same as next(int, str, str) but doesn't count sequence length.
                   Counting sequence length could be a very small overhead -> no need for this overload?

   @param  inext    IN     index of the stream to read.
   @param  readstr  OUT    read sequence
   @return true if record exists, else false
 */
bool Readfeed::next(int inext, std::string& readstr, bool is_orig)
{
	std::string line;
	std::stringstream readss; // an empty read
	auto stat = vstate_in[inext].last_stat;

	// read lines from the reads file and extract a single read
	for (auto count = vstate_in[inext].last_count; !vstate_in[inext].is_done; ++count) // count lines in a single record/read
	{
		if (vstate_in[inext].last_header.size() > 0)
		{
			if (!is_orig)
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
				if (!is_orig)
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
				if (count == 2) {
					if (is_orig)
						readss << line << '\n'; // line[0] == '+'
					continue;
				}
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
}

void Readfeed::rewind() {
	rewind_in();
}

/*
  rewind IN feed
*/
void Readfeed::rewind_in() {
	if (ifsv.size() >= num_orig_files) {
		for (auto i = 0; i < ifsv.size(); ++i) {
			if (ifsv[i].is_open()) {
				if (ifsv[i].rdstate() != std::ios_base::goodbit) {
					ifsv[i].clear();
				}
				ifsv[i].seekg(0); // rewind

				if (!ifsv[i].good()) {
					ERR("failed rewind stream: ", readfiles[i], " iostate: ", ifsv[i].rdstate());
					exit(1);
				}
			}
			vstate_in[i].reset();
		}
	}
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
	if (is_ready) {
		INFO("split is ready - no need to run");
		return true;
	}

	auto starts = std::chrono::high_resolution_clock::now();
	INFO("start splitting. Using number of splits equals number of processing threads: ", num_parts);
	auto retval = true;
	num_splits = num_parts;
	auto num_split_files = num_splits * num_orig_files;
	auto tot_files = num_split_files + num_orig_files; // splits + original read files

	// orig streams are ready at this stage - just check them
	if (ifsv.size() >= num_orig_files) {
		for (int i = 0; i < num_orig_files; ++i) {
			if (!ifsv[i].is_open()) {
				ifsv[i].open(readfiles[i], std::ios_base::in | std::ios_base::binary);
			}
			if (!ifsv[i].is_open()) {
				ERR("Failed to open file ", readfiles[i]);
				exit(1);
			}
		}
	}
	rewind_in();

	// init zlib for inflation
	for (int i = 0; i < num_orig_files; ++i) {
		if (vstate_in[i].isZip) {
			vzlib_in[i].init();
		}
	}

	// prepare split files to output:
	// [ fwd_1.fq.gz, rev_1.fq.gz,   fwd_2.fq.gz, rev_2.fq.gz,   ...,   fwd_n.fq.gz, ref_n.fq.gz ]
	ofsv.resize(num_split_files);
	readfiles.reserve(tot_files); // reserve space for new (split) names
	size_t idx = 0; // stream index
	std::stringstream ss;
	for (int i = 0; i < num_splits; ++i) {
		for (int j = 0; j < num_orig_files; ++j) {
			std::string stem = j == 0 ? "fwd_" : "rev_";
			std::string sfx_1 = vstate_in[j].isFasta ? ".fa" : ".fq";
			auto sfx = vstate_in[j].isZip ? sfx_1 + ".gz" : sfx_1;
			ss << stem << i << sfx; // split file basename
			auto fn = basedir / ss.str(); // split file name
			INFO("adding file: ", fn.generic_string());
			readfiles.emplace_back(fn.generic_string()); // add split file
			ss.str("");
		}
	}

	for (auto i = 0; i < ofsv.size(); ++i) {
		if (!ofsv[i].is_open()) {
			ofsv[i].open(readfiles[(size_t)i + num_orig_files], std::ios::app | std::ios::binary);
		}
		if (!ofsv[i].is_open()) {
			ERR("Failed to open file ", readfiles[(size_t)i + num_orig_files]);
			exit(1);
		}
		//INFO("opened file: ", readfiles[i + num_orig_files]);
	}

	// prepare zlib interface for writing split files
	vzlib_out.resize(num_split_files, Izlib(true, true));
	for (auto i = 0; i < vzlib_out.size(); ++i) {
		vzlib_out[i].init(true);
	}

	// prepare Readstates OUT
	vstate_out.resize(num_split_files);
	for (auto i = 0; i < num_splits; ++i) {
		auto del = i * num_orig_files;
		for (auto j = 0; j < num_orig_files; ++j) {
			auto idx = j + del;
			vstate_out[idx].isFasta = vstate_in[j].isFasta;
			vstate_out[idx].isFastq = vstate_in[j].isFastq;
			vstate_out[idx].isZip = vstate_in[j].isZip;
		}
	}

	// calculate number of reads in each of the output files
	auto nreads = num_orig_files == 2 ? num_reads_tot / 2 : num_reads_tot; // num reads in a single input file e.g. FWD
	auto minr = nreads / num_parts; // quotient i.e. min number of reads in each output file
	auto surplus = nreads - minr * num_parts; // remainder of reads to be distributed between the output files
	for (auto i = 0; i < num_parts; ++i) {
		auto maxr = i < surplus ? minr + 1 : minr; // distribute the surplus
		for (auto j = 0; j < num_orig_files; ++j) {
			vstate_out[j+i*num_orig_files].max_reads = maxr;
		}
	}

	unsigned seqlen = 0;
	std::string readstr;

	// loop until EOF - get reads - write into split files
	for (int inext = 0, isplit = 0, iout = 0, del = 0; next(inext, seqlen, readstr, true);)
	{
		++vstate_out[iout].read_count;
		auto is_last = vstate_out[iout].read_count == vstate_out[iout].max_reads;
		auto ret = vzlib_out[iout].defstr(readstr, ofsv[iout], is_last); // Z_STREAM_END | Z_OK - ok
		if (ret < Z_OK || ret > Z_STREAM_END) {
			ERR("Failed deflating readstring: ", readstr, " Output file idx: ", iout, " zlib status: ", ret);
			retval = false;
			break;
		}

		// split fwd, split rev:  FDW -> FDW_1,FDW_2,FWD_n; REV -> REV_1,REV_2,REV_n;
		if (is_last && inext == 1) {
			++isplit; // next split index
			del = isplit * num_orig_files; // distance between iout of 0th split and a current split
		}
		inext = is_two_files ? inext ^ 1 : inext; // toggle the index of the input file
		iout =  inext + del; // next out index
	} // ~for

	// close IN file streams
	for (int i = 0; i < num_orig_files; ++i) {
		if (ifsv[i].is_open()) {
			ifsv[i].close();
		}
	}

	// close Out streams
	for (auto i = 0; i < ofsv.size(); ++i) {
		if (ofsv[i].is_open())
			ofsv[i].close();
	}

	// reset Readstates. No need for izlib - already cleaned
	for (int i = 0; i < num_orig_files; ++i) {
		vstate_in[i].reset();
	}

	// zlib deflate streams already cleaned by now

	// write readfead descriptor

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

		// reset Izlib
		if (vstate_in[i].isZip) {
			// seems this Has to be done here i.e. within 
			// the same context where 'init' was called.
			// Otherwise 'Z_STREAM_ERROR'
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

	rewind_in();

	// init zlib
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
		if (max_read_len < seqlen) max_read_len = seqlen;
		if (min_read_len > seqlen || min_read_len == 0) min_read_len = seqlen;
		inext = is_two_files ? inext ^ 1 : inext; // toggle the index of the input file
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("done count. Elapsed time: ", elapsed.count(), " sec. Total reads: ", num_reads_tot);

} // ~Readfeed::count_reads

/*
  resize input streams vector to fit all split streams and open the streams
  resize and init Readstats and Izlibs
*/
void Readfeed::init_split_input() 
{
}

/*
  time:
  num_splits: 3
  num_orig_files: 2
  file:
  lines:
  file:
  lines:
  files:
  lines:
*/
void Readfeed::write_descriptor()
{
	std::filesystem::path fn = basedir / "readfeed";

}

void Readfeed::read_descriptor()
{
	auto nfiles = num_splits * num_orig_files;
	ifsv.resize((nfiles));
	auto ridx = num_orig_files; // split files follow original files in the 'ifsv'
	for (auto i = 0; i < ifsv.size(); ++i, ++ridx) {
		ifsv[i].open(readfiles[ridx], std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", readfiles[ridx]);
			exit(1);
		}
	}

	// prepare split Readstates
	vstate_in.resize(nfiles);
	vstate_in.resize(nfiles);
	for (auto i = 0; i < nfiles; ++i) {
		vstate_in[i].isFasta = vstate_out[i].isFasta;
		vstate_in[i].isFastq = vstate_out[i].isFastq;
		vstate_in[i].isZip = true;
		vstate_in[i].max_reads = vstate_out[i].max_reads;
	}

	// prepare Izlib
	vzlib_in.resize(nfiles);
	for (size_t i = 0; i < nfiles; ++i) {
		vzlib_in[i].init();
	}
} // ~Readfeed::read_descriptor

void Readfeed::init_vzlib_in()
{
	for (auto i = 0; i < vzlib_in.size(); ++i) {
		if (vstate_in[i].isZip)
			vzlib_in[i].init();
	}
}