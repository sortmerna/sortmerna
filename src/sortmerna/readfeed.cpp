/**
 * FILE: reader.cpp
 * Created: Nov 26, 2017 Sun
 * @copyright 2016-20 Clarity Genomics BVBA
 * 
 * Processes Reads file, creates Read objects, and pushes them to a queue for further pick−up by Processor
 */

#include <vector>
#include <iostream>
#include <ios> // std::ios_base
#include <algorithm> // find, find_if
#include <chrono> // std::chrono
#include <iomanip> // std::precision
#include <locale> // std::isspace
#include <thread>

#include "readfeed.hpp"

// forward
std::streampos filesize(const std::string& file); //util.cpp

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
	basedir(basedir)
{
	init(readfiles);
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
	basedir(basedir)
{
	init(readfiles);
} //~Readfeed::Readfeed 2

//Readfeed::~Readfeed() {}

void Readfeed::init(std::vector<std::string>& readfiles)
{
	// init read files
	orig_files.resize(num_orig_files);
	for (auto i = 0; i < num_orig_files; ++i) {
		orig_files[i].path = readfiles[i];
		orig_files[i].size = filesize(readfiles[i]);
	}

	// always do even when split is ready (need for validation)
	define_format();
	count_reads();

	if (type == FEED_TYPE::SPLIT_READS) {
		init_split_files();
		is_ready = is_split_ready();
		if (is_ready) { INFO("split is ready - no need to run"); }
		else split();
	}
	else {
		is_ready = true;
	}
}

/* 
 * can be run in a thread
 */
void Readfeed::run()
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("Readsfile::run thread ", std::this_thread::get_id(), " started");

	// init input file streams
	ifsv.resize(split_files.size());
	for (int i = 0; i < split_files.size(); ++i) {
		ifsv[i].open(split_files[i].path, std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", split_files[i].path);
			exit(1);
		}
	}

	// bug 114 - Z_STREAM_ERROR even though init is called at construction time
	vzlib_in.resize(split_files.size(), Izlib(true));
	for (auto i = 0; i < vzlib_in.size(); ++i) {
		vzlib_in[i].init(false);
	}

	std::string readstr;
	int inext = 0;
	unsigned seqlen = 0;

	// loop until EOF - get reads - push on queue
	for (;;)
	{
		if (next(inext, readstr, seqlen, split_files)) {
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

 @param inext   IN   index into the Readfeed::files vector
 @param seqlen  OUT  length of the read sequence
 @param read    OUT  read string
 @param is_orig IN   flags to return the original read string. If false, then the read string has format: 'read_id \n header \n sequence [\n quality]'
*/
bool Readfeed::next(int inext, std::string& read, unsigned& seqlen, bool is_orig, std::vector<Readfile>& files) {
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
			if (files[inext].isZip) {
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
			auto FR = (inext & 1) == 0 ? FWD : REV;
			INFO("EOF ", FR, " reached. Total reads: ", ++vstate_in[inext].read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			auto FR = (inext & 1) == 0 ? FWD : REV;
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
			files[inext].isFastq = (line[0] == FASTQ_HEADER_START);
			files[inext].isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && files[inext].isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((files[inext].isFasta && line[0] == FASTA_HEADER_START) || (files[inext].isFastq && count == 0)) // header line reached
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
			if (files[inext].isFastq)
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
   20201012: TODO: exactly the same as next(int, uint, str, bool) but doesn't count sequence length.
                   Counting sequence length could be a very small overhead -> no need for this overload?

   @param  inext    IN     index of the stream to read.
   @param  readstr  OUT    read sequence
   @return true if record exists, else false
 */
bool Readfeed::next(int inext, std::string& readstr, bool is_orig, std::vector<Readfile>& files)
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
			if (files[inext].isZip) {
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
			auto FR = (inext & 1) == 0 ? FWD : REV; // if index is even -> FWD
			INFO("EOF ", FR, " reached. Total reads: ", ++vstate_in[inext].read_count);
			break;
		}

		if (stat == RL_ERR)
		{
			auto FR = (inext & 1) == 0 ? FWD : REV;
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
			files[inext].isFastq = (line[0] == FASTQ_HEADER_START);
			files[inext].isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && files[inext].isFastq)
		{
			count = 0;
		}

		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if ((files[inext].isFasta && line[0] == FASTA_HEADER_START) || (files[inext].isFastq && count == 0)) // header line reached
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
			if (files[inext].isFastq)
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

	auto is_read_ok = readss.str().size() > 0;
	if (is_read_ok)
		readstr = readss.str();
	return is_read_ok;
} // ~Readfeed::next

/*
 * public function
 */
bool Readfeed::next(int inext, std::string& readstr)
{
	auto has_read = false;
	if (type == FEED_TYPE::SPLIT_READS)
		has_read = next(inext, readstr, false, split_files);
	return has_read;
}

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
	for (auto i = 0; i < ifsv.size(); ++i) {
		if (ifsv[i].is_open()) {
			if (ifsv[i].rdstate() != std::ios_base::goodbit) {
				ifsv[i].clear();
			}
			ifsv[i].seekg(0); // rewind

			if (!ifsv[i].good()) {
				ERR("failed rewind stream idx: ", i, " in vector of size: ", ifsv.size(), " iostate: ", ifsv[i].rdstate());
				exit(1);
			}
		}
		vstate_in[i].reset();
	}
}

/*
 input options:
   2 paired file     -> num files out = 2 x num_parts i.e. each paired file is split into specified number of parts
   1 paired file     -> num files out = 1 x num_parts AND ensure each part has even number of reads
   1 non-paired file -> num files out = 1 x num_parts
*/
bool Readfeed::split()
{
	if (is_ready) {
		INFO("split is ready - no need to run");
		return true;
	}

	auto starts = std::chrono::high_resolution_clock::now();
	INFO("start splitting. Using number of splits equals number of processing threads: ", num_splits);
	auto retval = true;
	auto num_split_files = num_splits * num_orig_files;

	// orig streams are ready at this stage - just check them
	for (auto i = 0; i < ifsv.size(); ++i) {
		if (!ifsv[i].is_open()) {
			ifsv[i].open(orig_files[i].path, std::ios_base::in | std::ios_base::binary);
		}
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", orig_files[i].path);
			exit(1);
		}
	}
	rewind_in();

	// init zlib IN for inflation
	for (int i = 0; i < orig_files.size(); ++i) {
		if (orig_files[i].isZip) {
			vzlib_in[i].init();
		}
	}

	// prepare split files to output:
	//init_split_files(num_splits);
	// [ fwd_1.fq.gz, rev_1.fq.gz,   fwd_2.fq.gz, rev_2.fq.gz,   ...,   fwd_n.fq.gz, ref_n.fq.gz ]
	ofsv.resize(num_split_files);

	for (auto i = 0; i < ofsv.size(); ++i) {
		if (!ofsv[i].is_open()) {
			ofsv[i].open(split_files[i].path, std::ios::app | std::ios::binary);
		}
		if (!ofsv[i].is_open()) {
			auto idx = num_orig_files + i;
			ERR("Failed to open file ", split_files[i].path.generic_string());
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

	unsigned seqlen = 0;
	std::string readstr;

	// loop until EOF - get reads - write into split files
	for (auto inext = 0, iout = 0; inext < num_orig_files; ++inext, iout = inext) {
		for (; next(inext, readstr, seqlen, true, orig_files);)
		{
			++vstate_out[iout].read_count;
			auto ret = vzlib_out[iout].defstr(readstr, ofsv[iout], vstate_out[iout].read_count == split_files[iout].numreads); // Z_STREAM_END | Z_OK - ok
			if (ret < Z_OK || ret > Z_STREAM_END) {
				ERR("Failed deflating readstring: ", readstr, " Output file idx: ", iout, " zlib status: ", ret);
				retval = false;
				break;
			}
			if (vstate_out[iout].read_count == split_files[iout].numreads) {
				iout += num_orig_files; // switch to next split
			}
		} // ~for
	}

	// loop until EOF - get reads - write into split files
	//for (int inext = 0, isplit = 0, iout = 0, del = 0; next(inext, readstr, seqlen, true, orig_files);)
	//{
	//	++vstate_out[iout].read_count;
	//	auto ofdx = iout + num_orig_files; // split file index into files vector
	//	auto is_last = vstate_out[iout].read_count == split_files[ofdx].numreads;
	//	auto ret = vzlib_out[iout].defstr(readstr, ofsv[iout], is_last); // Z_STREAM_END | Z_OK - ok
	//	if (ret < Z_OK || ret > Z_STREAM_END) {
	//		ERR("Failed deflating readstring: ", readstr, " Output file idx: ", iout, " zlib status: ", ret);
	//		retval = false;
	//		break;
	//	}
	//
	//	if (is_last && inext == 1) {
	//		++isplit; // next split index
	//		del = isplit * num_orig_files; // distance between iout of 0th split and a current split
	//	}
	//	inext = is_two_files ? inext ^ 1 : inext; // toggle the index of the input file
	//	iout =  inext + del; // next out index
	//} // ~for

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

	// calculate split file sizes
	for (auto i = 0; i < split_files.size(); ++i) {
		split_files[i].size = filesize(split_files[i].path.generic_string());
	}

	// write readfead descriptor
	write_descriptor();

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("Done splitting. Reads count: ", num_reads_tot, " Runtime sec: ", elapsed.count(), "\n");
	return retval;
} // ~Readfeed::split

/*
  verify the split was already performed and the feed is ready
  logic:
	Read feed descriptor and verify:
	  - original files have the same name, size and line count
	  - num splits are the same
	  - split files names, count, and line count are the same
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
	// compare with data in descriptor
	std::ifstream ifs;
	std::filesystem::path fn = basedir / "readfeed";
	if (std::filesystem::exists(fn)) {
		ifs.open(fn, std::ios_base::in | std::ios_base::binary);
		if (!ifs.is_open()) {
			ERR("failed to open: ", fn.generic_string());
			exit(1);
		}

		int lidx = 0; // line index
		int fcnt = 0; // count of file entries in the desctiptor
		int fcnt_max = (1 + num_splits) * num_orig_files;
		int fidx = 0; // file index
		int fpidx = 0; // index of file parameters: name, size, lines, zip, fastq/fasta
		for (std::string line; std::getline(ifs, line); ) {
			if (!line.empty()) {
				// trim
				line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
				if (line[0] != '#') { // skip comments
					if (lidx == 0)
						is_ready = true; // timestamp
					else if (lidx == 1)
						is_ready = is_ready && std::stoi(line) == num_orig_files;
					else if (lidx == 2)
						is_ready = is_ready && std::stoi(line) == num_splits;
					else if (lidx == 3)
						is_ready = is_ready && std::stoi(line) == num_reads_tot;
					else {
						if (fpidx == 0) { // file path
							++fcnt;
							if (fcnt == num_orig_files+1)	fidx = 0; // switch file index orig -> split files
							if (fcnt > fcnt_max) {
								INFO("expected max file count: ", fcnt_max, " current file count in the descriptor: ", fcnt);
								is_ready = false;
								break;
							}
							if (fcnt <= num_orig_files)
								is_ready = is_ready && line == orig_files[fidx].path;
							else
								is_ready = is_ready && line == split_files[fidx].path;
						}
						else if (fpidx == 1) { // file size
							if (fcnt <= num_orig_files)
								is_ready = is_ready && std::stoi(line) == orig_files[fidx].size; // original files
							else
								split_files[fidx].size = std::stoi(line); // set sizes of split files
						}
						else if (fpidx == 2) { // number of reads
							if (fcnt <= num_orig_files)
								is_ready = is_ready && std::stoi(line) == orig_files[fidx].numreads;
							else
								is_ready = is_ready && std::stoi(line) == split_files[fidx].numreads;
						}
						else if (fpidx == 3) { // is zip
							auto isZip = std::stoi(line) == 1;
							if (fcnt <= num_orig_files)
								is_ready = is_ready && orig_files[fidx].isZip == isZip;
							else
								is_ready = is_ready && split_files[fidx].isZip == isZip;
						}
						else if (fpidx == 4) { // fastq/a
							auto isFastq = line == "fastq";
							auto isFasta = line == "fasta";
							if (fcnt <= num_orig_files) {
								is_ready = is_ready && orig_files[fidx].isFastq == isFastq;
								is_ready = is_ready && orig_files[fidx].isFasta == isFasta;
							}
							else {
								is_ready = is_ready && split_files[fidx].isFastq == isFastq;
								is_ready = is_ready && split_files[fidx].isFasta == isFasta;
							}

							fpidx = 0;
							++fidx;
							++lidx;
							continue;
						}
						++fpidx;
					}
					if (!is_ready)
						break;
					++lidx;
				}
			}
		} // ~for lines in descriptor
		// verify count of files
		is_ready = is_ready && (size_t)fidx == split_files.size();

		// reset split file sizes if not ready
		if (!is_ready) {
			for (auto i = 0; i < split_files.size(); ++i) {
				split_files[i].size = 0;
			}
		}
	}
	else {
		is_ready = false; // descriptor does not exist
	}

	if (ifs.is_open()) ifs.close(); // close the descriptor

	return is_ready;
} // ~Readfeed::is_split_ready

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

	ifsv.resize(num_orig_files);
	vstate_in.resize(num_orig_files);

	std::string str(100, '\0');
	unsigned seqlen = 0;
	for (int i = 0; i < num_orig_files; ++i) {
		if (!ifsv[i].is_open()) {
			ifsv[i].open(orig_files[i].path, std::ios_base::in | std::ios_base::binary);
		}
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", orig_files[i].path);
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
		orig_files[i].isZip = !is_ascii;
		if (!is_ascii) {
			// init izlib for inflation
			if (vzlib_in.size() < i+1) {
				vzlib_in.emplace_back(Izlib());
			}
			vzlib_in[i].init(false);
		}
		
		ifsv[i].seekg(0); // rewind stream to the start
		// get a read
		if (next(i, str, seqlen, orig_files)) {
			std::string fmt = orig_files[i].isZip ? "gzipped" : "flat ASCII";
			if (orig_files[i].isFasta) {
				biof = BIO_FORMAT::FASTA;
				INFO("file: ", orig_files[i].path, " is FASTA ", fmt);
			}
			else if (orig_files[i].isFastq) {
				biof = BIO_FORMAT::FASTQ;
				INFO("file: ", orig_files[i].path, " is FASTQ ", fmt);
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
			ERR("Cannot define format for file: ", orig_files[i].path);
			exit(1);
		}

		// reset Izlib
		if (orig_files[i].isZip) {
			// seems this Has to be done here i.e. within 
			// the same context where 'init' was called.
			// Otherwise 'Z_STREAM_ERROR'
			auto stat = vzlib_in[i].reset_inflate();
		}
	} // ~for orig files

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
	for (int i = 0; i < orig_files.size(); ++i) {
		if (orig_files[i].isZip) {
			vzlib_in[i].init();
		}
	}

	std::string readstr;
	unsigned seqlen = 0;

	// loop until EOF - count reads and sequence lengths
	for (int inext = 0; next(inext, readstr, seqlen, true, orig_files);)
	{
		++num_reads_tot;
		++orig_files[inext].numreads;
		length_all += seqlen;
		if (max_read_len < seqlen) max_read_len = seqlen;
		if (min_read_len > seqlen || min_read_len == 0) min_read_len = seqlen;
		inext = is_two_files ? inext ^ 1 : inext; // toggle the index of the input file
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("done count. Elapsed time: ", elapsed.count(), " sec. Total reads: ", num_reads_tot);

} // ~Readfeed::count_reads

/*
  init split files' names, and attributes including format fastq/a, zip, number of reads.
*/
void Readfeed::init_split_files()
{
	auto n_split_files = num_orig_files * num_splits;
	split_files.resize(n_split_files);

	std::stringstream ss;
	for (int i = 0, idx = 0; i < num_splits; ++i) {
		for (int j = 0; j < num_orig_files; ++j, ++idx) {
			std::string stem = j == 0 ? "fwd_" : "rev_";
			std::string sfx_1 = orig_files[j].isFasta ? ".fa" : ".fq";
			auto sfx = orig_files[j].isZip ? sfx_1 + ".gz" : sfx_1;
			ss << stem << i << sfx; // split file basename
			//auto idx = i * num_orig_files + j;
			split_files[idx].path = basedir / ss.str();
			ss.str("");
			INFO("added file: ", split_files[idx].path.generic_string());
		}
	}

	// calculate number of reads in each of the split files
	auto nreads = num_orig_files == 2 ? num_reads_tot / 2 : num_reads_tot; // num reads in a single input file e.g. FWD
	auto minr = nreads / num_splits; // quotient i.e. min number of reads in each output file
	auto surplus = nreads - minr * num_splits; // remainder of reads to be distributed between the output files
	for (auto i = 0, idx = 0; i < num_splits; ++i) {
		auto maxr = i < surplus ? minr + 1 : minr; // distribute the surplus
		for (auto j = 0; j < num_orig_files; ++j, ++idx) {
			//auto idx = i * num_orig_files + num_orig_files + j;
			split_files[idx].numreads = maxr;
			split_files[idx].isZip = orig_files[j].isZip;
			split_files[idx].isFastq = orig_files[j].isFastq;
			split_files[idx].isFasta = orig_files[j].isFasta;
		}
	}
}

/*
  format:
    time:
    num_orig_files: 2
    num_splits: 3
	total_reads: xx
    [ file:
      size:
      reads:
    ]  <- as many as num_orig_files*(num_splits + 1)
*/
void Readfeed::write_descriptor()
{
	const std::string comments = 
		"# format of this file:\n"
	    "#   time\n"
		"#   num_orig_files\n"
		"#   num_splits\n"
		"#   num_reads_tot\n"
        "#   [\n"
		"#     file\n"
		"#     size\n"
		"#     reads\n"
		"#     zip\n"
		"#     fastq/a\n"
        "#   ] for each file both original and split\n";

	if (!is_ready) {
		std::filesystem::path fn = basedir / "readfeed";
		std::ofstream ofs(fn, std::ios::app | std::ios::binary);
		if (!ofs.is_open()) {
			ERR("failed to open file: ", fn.generic_string());
			exit(1);
		}

		ofs << comments << '\n';
		std::time_t tm = std::time(0);
		ofs << std::ctime(&tm) << '\n'; // line 0: timestamp  <ctime> 'Tue Oct 20 08:39:35 2020'
		ofs << num_orig_files << '\n';  // line 1:
		ofs << num_splits << '\n';      // line 2:
		ofs << num_reads_tot << '\n';   // line 3:
		for (auto i = 0; i < orig_files.size(); ++i) {
			ofs << orig_files[i].path.generic_string() << '\n';     // file path
			ofs << orig_files[i].size << '\n';     // file size
			ofs << orig_files[i].numreads << '\n'; // reads count
			ofs << orig_files[i].isZip << '\n';    // zip
			auto fx = orig_files[i].isFastq == true ? "fastq" : "fasta";
			ofs << fx << '\n';
		}
		for (auto i = 0; i < split_files.size(); ++i) {
			ofs << split_files[i].path.generic_string() << '\n';     // file path
			ofs << split_files[i].size << '\n';     // file size
			ofs << split_files[i].numreads << '\n'; // reads count
			ofs << split_files[i].isZip << '\n';    // zip
			auto fx = split_files[i].isFastq == true ? "fastq" : "fasta";
			ofs << fx << '\n';
		}
		is_ready = true;
	}
	else {
		INFO("split is ready - no write_descriptor necessary");
	}
} // ~Readfeed::write_descriptor

void Readfeed::init_vzlib_in()
{
	vzlib_in.resize(split_files.size());
	for (auto i = 0; i < vzlib_in.size(); ++i) {
		if (split_files[i].isZip) vzlib_in[i].init();
	}
}

/*
 * called at the start of reading the split files 
 * init readfeed for reading:
 *   ifsv
 *   vstate_in
 *   vzlib_in
 */
void Readfeed::init_reading()
{
	for (auto i = 0; i < vstate_in.size(); ++i) {
		vstate_in[i].reset();
	}
	vstate_in.resize(split_files.size());

	// ifsv
	for (auto i = 0; i < ifsv.size(); ++i) {
		if (ifsv[i].is_open()) ifsv[i].close();
	}
	ifsv.resize(split_files.size());
	for (auto i = 0; i < split_files.size(); ++i) {
		ifsv[i].open(split_files[i].path, std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("failed to open: ", split_files[i].path.generic_string());
			exit(1);
		}
	}

	// vzlib_in
	init_vzlib_in();
} // ~Readfeed::init_reading