/*
@copyright 2016-2026 Clarity Genomics BVBA
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
              biocodz          biocodz@protonmail.com
*/

/*
 * FILE: readfeed.cpp
 * Created: Nov 26, 2017 Sun
 * 
 */

#include "common.hpp"
#include "readfeed.hpp"

#include <vector>
#include <iostream>
#include <ios> // std::ios_base
#include <algorithm> // find, find_if
#include <chrono> // std::chrono
#include <iomanip> // std::precision
#include <locale> // std::isspace
#include <thread>
#include <regex>

#include <filereader/Standard.hpp>
#include <rapidgzip/ParallelGzipReader.hpp>

// Opaque wrapper — keeps rapidgzip headers out of readfeed.hpp and every TU that includes it.
// The explicit template specialisation NEXT_DYNAMIC_DEFLATE_CANDIDATE_LUT<15> in DynamicHuffman.hpp
// lacks 'inline', so it gets external linkage; Apple ld rejects duplicate definitions across TUs.
struct GzReaderImpl {
    rapidgzip::ParallelGzipReader<> rdr;
    template <typename T>
    GzReaderImpl(std::unique_ptr<T> fr, std::size_t threads)
        : rdr(std::move(fr), threads) {}
};

// Custom deleter — body defined here so sizeof(GzReaderImpl) is never checked in other TUs.
void GzReaderDeleter::operator()(GzReaderImpl* p) noexcept { delete p; }

// forward
std::streampos filesize(const std::string& file); //util.cpp

// ---------------------------------------------------------------------------
// FlatSlot implementation
// ---------------------------------------------------------------------------

bool FlatSlot::fill_buf()
{
	if (bytes_remaining == 0) return false;
	size_t toRead = static_cast<size_t>(std::min(static_cast<uint64_t>(BUF_SIZE), bytes_remaining));
	ifs.read(buf.data(), static_cast<std::streamsize>(toRead));
	auto n = ifs.gcount();
	if (n <= 0) { bytes_remaining = 0; return false; }
	bytes_remaining -= static_cast<uint64_t>(n);
	buf_pos = 0;
	buf_len = static_cast<size_t>(n);
	return true;
}

int FlatSlot::getline(std::string& line)
{
	line.clear();
	for (;;) {
		if (buf_pos >= buf_len) {
			if (!fill_buf())
				return line.empty() ? RL_END : RL_OK;
		}
		while (buf_pos < buf_len) {
			char c = buf[buf_pos++];
			if (c == '\n') return RL_OK;
			line += c;
		}
	}
}

// ---------------------------------------------------------------------------
// GzSlot implementation
// ---------------------------------------------------------------------------

bool GzSlot::fill_buf()
{
	if (bytes_remaining == 0) return false;
	size_t toRead = static_cast<size_t>(std::min(static_cast<uint64_t>(BUF_SIZE), bytes_remaining));
	auto n = reader->rdr.read(reinterpret_cast<char*>(buf.data()), toRead);
	if (n <= 0) { bytes_remaining = 0; return false; }
	bytes_remaining -= static_cast<uint64_t>(n);
	buf_pos = 0;
	buf_len = static_cast<size_t>(n);
	return true;
}

int GzSlot::getline(std::string& line)
{
	line.clear();
	for (;;) {
		if (buf_pos >= buf_len) {
			if (!fill_buf())
				return line.empty() ? RL_END : RL_OK; // EOF or end-of-chunk
		}
		while (buf_pos < buf_len) {
			char c = static_cast<char>(buf[buf_pos++]);
			if (c == '\n') return RL_OK;
			line += c;
		}
	}
}

/*
 @param type       feed type
 @param readfiles  vector with reads file paths
 @param is_gz      flags the readsfiles format gzip | non-gzip (flat)

 Split reads logic:
   - check the split was already done:
     - check readb directory for the descriptor and file
*/
Readfeed::Readfeed(FEED_TYPE type, 
                    std::vector<std::string>& readfiles, 
                    std::filesystem::path& basedir, 
                    bool is_paired)
	:
	type(type),
	is_done(false),
	is_ready(false),
	is_format_defined(false),
	is_two_files(readfiles.size() > 1),
	is_paired(is_paired),
	num_orig_files(readfiles.size()),
	num_splits(0),
	num_split_files(0),
	num_sense(0),
	num_reads_tot(0),
	length_all(0),
	min_read_len(0),
	max_read_len(0),
	basedir(basedir)
{
	init(readfiles);
} // ~Readfeed::Readfeed 1

Readfeed::Readfeed(FEED_TYPE type, 
                    std::vector<std::string>& readfiles, 
                    const unsigned num_parts, 
                    std::filesystem::path& basedir, 
                    bool is_paired)
	:
	type(type),
	is_done(false),
	is_ready(false),
	is_format_defined(false),
	is_two_files(readfiles.size() > 1),
	is_paired(is_paired),
	num_orig_files(readfiles.size()),
	num_splits(num_parts),
	num_split_files(0),
	num_sense(0),
	num_reads_tot(0),
	length_all(0),
	min_read_len(0),
	max_read_len(0),
	basedir(basedir)
{
	init(readfiles);
} //~Readfeed::Readfeed 2

//Readfeed::~Readfeed() {}

void Readfeed::init(std::vector<std::string>& readfiles, const int& dbg)
{
	auto start = std::chrono::high_resolution_clock::now();
	INFO("Readfeed init started");

	num_sense = is_paired ? 2 : 1;
	num_split_files = num_sense * num_splits;

	// init read files
	orig_files.resize(num_orig_files);
	for (decltype(num_orig_files) i = 0; i < num_orig_files; ++i) {
		orig_files[i].path = readfiles[i];
		orig_files[i].size = filesize(readfiles[i]);
	}

	define_format();
    // calculate this.num_reads_tot
    if (type == FEED_TYPE::INDEXED) {
        count_reads_parallel();
    }
    else {
        count_reads();
    }

	// Select INDEXED_GZ or INDEXED_FLAT based on input file format.
    // Single interleaved paired files when (num_orig_files < num_sense) 
	// are handled by sharing a single reader between FWD/REV slots.
	//if (num_splits > 0) {
	//	bool all_gz_fastq = true;
	//	bool all_flat = true;
	//	for (auto& f : orig_files) {
	//		if (!f.isZip || !f.isFastq) all_gz_fastq = false;
	//		if (f.isZip)                all_flat     = false;
	//	}
	//	const bool is_interleaved = (num_orig_files < num_sense);
	//	if (all_gz_fastq) {
	//		if (is_interleaved) { 
    //            INFO("Input: single interleaved gzipped FASTQ — using INDEXED_GZ feed type"); 
    //        }
	//		else { 
    //            INFO("Input: gzipped FASTQ — using INDEXED_GZ feed type"); 
    //        }
	//		type = FEED_TYPE::INDEXED_GZ;
	//	} else if (all_flat) {
	//		if (is_interleaved) { INFO("Input: single interleaved flat — using INDEXED_FLAT feed type"); }
	//		else                { INFO("Input: flat — using INDEXED_FLAT feed type"); }
	//		type = FEED_TYPE::INDEXED_FLAT;
	//	} else {
	//		WARN("Mixed gzipped/flat input — falling back to deprecated SPLIT_READS feed type");
	//	}
	//}
	if (type == FEED_TYPE::INDEXED) {
        if (orig_files[0].isZip) {
		    build_chunk_offsets();
        }
        else {
		    build_flat_chunk_offsets();
        }
		is_ready = true;
	}
	else if (type == FEED_TYPE::SPLIT_READS) {
		init_split_files();
		is_ready = is_split_ready();
		if (is_ready) { INFO("split is ready - no need to run"); }
		else split();
	}
	else {
        // should never get here since feed type is validated at the command line parsing stage, but just in case...
        ERR("Unsupported feed type: ", static_cast<unsigned>(type));
        exit(1);
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("Readfeed init done in sec [", elapsed.count(), "]\n");
}

/* 
 * can be run in a threadq
 */
void Readfeed::run()
{
	auto starts = std::chrono::high_resolution_clock::now();
	INFO("Readsfile::run thread ", std::this_thread::get_id(), " started");

	// init input file streams
	ifsv.resize(split_files.size());
	for (std::size_t i = 0; i < split_files.size(); ++i) {
		ifsv[i].open(split_files[i].path, std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", split_files[i].path);
			exit(1);
		}
	}

	// bug 114 - Z_STREAM_ERROR even though init is called at construction time
	vzlib_in.resize(split_files.size(), Izlib(true));
	for (std::size_t i = 0; i < vzlib_in.size(); ++i) {
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
		for (std::size_t i = 0; i < vstate_in.size(); ++i) {
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

bool Readfeed::next(int inext, 
                    std::string& read, 
                    unsigned& seqlen, 
                    bool is_orig, 
                    std::vector<Readfile>& files) {
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
		line = "";
		if (!vstate_in[inext].is_done) {
			if (files[inext].isZip) {
				stat = vzlib_in[inext].getline(ifsv[inext], line); // zipped
			}
			else {
				if (ifsv[inext].eof())
					stat = RL_END;
				else 
					std::getline(ifsv[inext], line);  // non-zipped
			}
		}

		// right-trim whitespace in place (removes '\r' too)
		if (!line.empty()) {
			line.erase(std::find_if(line.rbegin(), line.rend(), 
                        [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
		}

		// EOF reached - return last read
		if (stat == RL_END)
		{
			// 20210714 - issue 294 - add last line (FQ quality or FA sequence)
			//         for cases where there is no NL at the end of reads file
			if (!line.empty()) {
				readss << line;
			}
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

		// define fasta/q on the first line in file
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
				// 20210714 - issue 294 - add NL to fasta. Cannot do earlier due possible multiline fasta sequence
				if (files[inext].isFasta)
					readss << '\n';
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
		line = "";
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

		if (!line.empty()) {
			// right-trim whitespace in place (removes '\r' too)
			line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
		}

		// EOF reached - return last read
		if (stat == RL_END)
		{
			// 20210714 - issue 294 - add last line (FQ quality or FA sequence)
			//         for cases where there is no NL at the end of reads file
			if (!line.empty()) {
				readss << line;
			}
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
				// 20210714 - issue 294 - add NL to fasta. Cannot do earlier due possible multiline fasta sequence
				if (files[inext].isFasta)
					readss << '\n';
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
 * next_gz  (INDEXED_GZ)
 * Mirrors next(int,string&,bool,vector<Readfile>&) but reads lines from
 * gz_slots[inext].getline() instead of vzlib_in / ifsv.
 */
bool Readfeed::next_gz(int inext, std::string& readstr, bool is_orig)
{
	// For interleaved paired (single file), FWD and REV slots share one reader.
	// REV slot (inext % num_sense != 0) delegates to its FWD partner's reader and state.
	const int slot_idx = (num_orig_files < num_sense && inext % static_cast<int>(num_sense) != 0)
	                     ? inext - (inext % static_cast<int>(num_sense))
	                     : inext;

	std::string line;
	std::stringstream readss;
	auto stat = vstate_in[slot_idx].last_stat;
	auto& files = gz_slot_files;

	for (auto count = vstate_in[slot_idx].last_count; !vstate_in[slot_idx].is_done; ++count)
	{
		if (vstate_in[slot_idx].last_header.size() > 0) {
			if (!is_orig)
				readss << inext << '_' << vstate_in[slot_idx].read_count << '\n';
			readss << vstate_in[slot_idx].last_header << '\n';
			vstate_in[slot_idx].last_header = "";
		}

		line = "";
		if (!vstate_in[slot_idx].is_done) {
			stat = gz_slots[slot_idx].getline(line);
		}

		if (!line.empty()) {
			line.erase(std::find_if(line.rbegin(), line.rend(),
				[l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
		}

		if (stat == RL_END) {
			if (!line.empty()) readss << line;
			vstate_in[slot_idx].is_done = true;
			auto FR = (inext & 1) == 0 ? FWD : REV;
			INFO("EOF ", FR, " reached. Total reads: ", ++vstate_in[slot_idx].read_count);
			break;
		}

		if (stat == RL_ERR) {
			auto FR = (inext & 1) == 0 ? FWD : REV;
			ERR("reading from ", FR, " file. Exiting...");
			exit(1);
		}

		if (line.empty()) { --count; continue; }

		++vstate_in[slot_idx].line_count;

		if (vstate_in[slot_idx].line_count == 1) {
			files[slot_idx].isFastq = (line[0] == FASTQ_HEADER_START);
			files[slot_idx].isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && files[slot_idx].isFastq) count = 0;

		if ((files[slot_idx].isFasta && line[0] == FASTA_HEADER_START) ||
			(files[slot_idx].isFastq && count == 0))
		{
			if (vstate_in[slot_idx].line_count == 1) {
				if (!is_orig)
					readss << inext << '_' << vstate_in[slot_idx].read_count << '\n';
				readss << line << '\n';
				count = 0;
			} else {
				if (files[slot_idx].isFasta) readss << '\n';
				vstate_in[slot_idx].last_header = line;
				vstate_in[slot_idx].last_count  = 1;
				vstate_in[slot_idx].last_stat   = stat;
				break;
			}
		} else {
			if (files[slot_idx].isFastq) {
				if (count == 2) { if (is_orig) readss << line << '\n'; continue; }
				if (count == 3) { readss << line; if (is_orig) readss << '\n'; continue; }
				readss << line << '\n';
			} else {
				readss << line;
			}
		}
	} // ~for getline

	++vstate_in[slot_idx].read_count;
	auto is_read_ok = readss.str().size() > 0;
	if (is_read_ok) readstr = readss.str();
	return is_read_ok;
} // ~Readfeed::next_gz

/*
 * next_flat  (INDEXED_FLAT)
 * Mirrors next_gz() but reads lines from flat_slots[inext].getline().
 */
bool Readfeed::next_flat(int inext, std::string& readstr, bool is_orig)
{
	// For interleaved paired (single file), FWD and REV slots share one ifstream.
	// REV slot (inext % num_sense != 0) delegates to its FWD partner's ifstream and state.
	const int slot_idx = (num_orig_files < num_sense && inext % static_cast<int>(num_sense) != 0)
	                     ? inext - (inext % static_cast<int>(num_sense))
	                     : inext;

	std::string line;
	std::stringstream readss;
	auto stat = vstate_in[slot_idx].last_stat;
	auto& files = flat_slot_files;

	for (auto count = vstate_in[slot_idx].last_count; !vstate_in[slot_idx].is_done; ++count)
	{
		if (vstate_in[slot_idx].last_header.size() > 0) {
			if (!is_orig)
				readss << slot_idx << '_' << vstate_in[slot_idx].read_count << '\n';
			readss << vstate_in[slot_idx].last_header << '\n';
			vstate_in[slot_idx].last_header = "";
		}

		line = "";
		if (!vstate_in[slot_idx].is_done) {
			stat = flat_slots[slot_idx].getline(line);
		}

		if (!line.empty()) {
			line.erase(std::find_if(line.rbegin(), line.rend(),
				[l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
		}

		if (stat == RL_END) {
			if (!line.empty()) readss << line;
			vstate_in[slot_idx].is_done = true;
            if (num_orig_files == 1) {
			    INFO("EOF reached. Slot: ", slot_idx, " Total reads: ", ++vstate_in[slot_idx].read_count);
            }
            else {
			    auto FR = (inext & 1) == 0 ? FWD : REV;
			    INFO("EOF ", FR, " reached. Slot: ", slot_idx, " Total reads: ", ++vstate_in[slot_idx].read_count);
            }
			break;
		}

		if (stat == RL_ERR) {
            if (num_orig_files == 1) {
                ERR("reading from file. Slot: ", slot_idx, " Exiting...");
            }
            else {
			    auto FR = (inext & 1) == 0 ? FWD : REV;
			    ERR("reading from ", FR, " file. Slot: ", slot_idx, " Exiting...");
            }
			exit(1);
		}

		if (line.empty()) { --count; continue; }

		++vstate_in[slot_idx].line_count;

		if (vstate_in[slot_idx].line_count == 1) {
			files[slot_idx].isFastq = (line[0] == FASTQ_HEADER_START);
			files[slot_idx].isFasta = (line[0] == FASTA_HEADER_START);
		}

		if (count == 4 && files[slot_idx].isFastq) count = 0;

		if ((files[slot_idx].isFasta && line[0] == FASTA_HEADER_START) ||
			(files[slot_idx].isFastq && count == 0))
		{
			if (vstate_in[slot_idx].line_count == 1) {
				if (!is_orig)
					readss << slot_idx << '_' << vstate_in[slot_idx].read_count << '\n';
				readss << line << '\n';
				count = 0;
			} else {
				if (files[slot_idx].isFasta) readss << '\n';
				vstate_in[slot_idx].last_header = line;
				vstate_in[slot_idx].last_count  = 1;
				vstate_in[slot_idx].last_stat   = stat;
				break;
			}
		} else {
			if (files[slot_idx].isFastq) {
				if (count == 2) { if (is_orig) readss << line << '\n'; continue; }
				if (count == 3) { readss << line; if (is_orig) readss << '\n'; continue; }
				readss << line << '\n';
			} else {
				readss << line;
			}
		}
	} // ~for getline

	++vstate_in[slot_idx].read_count;
	auto is_read_ok = readss.str().size() > 0;
	if (is_read_ok) readstr = readss.str();
	return is_read_ok;
} // ~Readfeed::next_flat

/*
 * public function
 */
bool Readfeed::next(int inext, std::string& readstr)
{
	if (type == FEED_TYPE::SPLIT_READS)
		return next(inext, readstr, false, split_files);
	if (type == FEED_TYPE::INDEXED && orig_files[0].isZip)
		return next_gz(inext, readstr, false);
	if (type == FEED_TYPE::INDEXED && !orig_files[0].isZip)
		return next_flat(inext, readstr, false);
	return false;
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
	if (type == FEED_TYPE::INDEXED && orig_files[0].isZip) {
		const bool is_interleaved = (num_orig_files < num_sense);
		for (std::size_t i = 0; i < gz_slots.size(); ++i) {
			if (is_interleaved && i % num_sense != 0) continue; // REV slots share FWD reader
			auto& slot = gz_slots[i];
			if (slot.reader) {
				slot.reader->rdr.seek(static_cast<long long>(slot.bytes_start));
			}
			slot.bytes_remaining = slot.bytes_end - slot.bytes_start;
			slot.buf_pos = 0;
			slot.buf_len = 0;
			if (i < vstate_in.size()) vstate_in[i].reset();
		}
		return;
	}

	if (type == FEED_TYPE::INDEXED && !orig_files[0].isZip) {
		const bool is_interleaved = (num_orig_files < num_sense);
		for (std::size_t i = 0; i < flat_slots.size(); ++i) {
			if (is_interleaved && i % num_sense != 0) continue; // REV slots share FWD ifstream
			auto& slot = flat_slots[i];
			if (slot.ifs.is_open()) {
				if (slot.ifs.rdstate() != std::ios_base::goodbit) slot.ifs.clear();
				slot.ifs.seekg(static_cast<std::streamoff>(slot.bytes_start));
			}
			slot.bytes_remaining = slot.bytes_end - slot.bytes_start;
			slot.buf_pos = 0;
			slot.buf_len = 0;
			if (i < vstate_in.size()) vstate_in[i].reset();
		}
		return;
	}

    // SPLIT_READS - deprecated
    // 20260329 Sun TODO: remove after INDEXED_GZ and INDEXED_FLAT are fully supported and tested
	for (std::size_t i = 0; i < ifsv.size(); ++i) {
		if (ifsv[i].is_open()) {
			if (ifsv[i].rdstate() != std::ios_base::goodbit) {
				ifsv[i].clear();
			}
			ifsv[i].seekg(0); // rewind

			if (!ifsv[i].good()) {
				ERR("failed rewind stream idx: ", i, " in vector of size: ", ifsv.size(), 
                    " iostate: ", ifsv[i].rdstate());
				exit(1);
			}
		}
		vstate_in[i].reset();
	}
}

/*
 * read input files (possibly compressed) and output split files (compressed)
 *
 * input options:
 *   2 paired file     -> num files out = 2 x num_parts i.e. each paired file is split into specified number of parts
 *   1 paired file     -> num files out = 1 x num_parts AND ensure each part has even number of reads
 *   1 non-paired file -> num files out = 1 x num_parts
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

	// remove existing split files
	clean();

	// orig streams are ready at this stage - just check them
	for (std::size_t i = 0; i < ifsv.size(); ++i) {
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
	bool is_zip = false; // flags at least one file requires archiving operations
	for (std::size_t i = 0; i < orig_files.size(); ++i) {
		if (orig_files[i].isZip) {
			vzlib_in[i].init();
			is_zip = true;
		}
	}

	// prepare split files to output:
	//init_split_files(num_splits);
	// [ fwd_1.fq.gz, rev_1.fq.gz,   fwd_2.fq.gz, rev_2.fq.gz,   ...,   fwd_n.fq.gz, ref_n.fq.gz ]
	ofsv.resize(num_split_files);

	for (std::size_t i = 0; i < ofsv.size(); ++i) {
		if (!ofsv[i].is_open()) {
			ofsv[i].open(split_files[i].path, std::ios::out | std::ios::binary | std::ios::trunc);
		}
		if (!ofsv[i].is_open()) {
			ERR("Failed to open file ", split_files[i].path.generic_string());
			exit(1);
		}
	}

	// prepare zlib interface for writing split files
	if (is_zip) {
		vzlib_out.resize(num_split_files, Izlib(true, true));
		for (std::size_t i = 0; i < vzlib_out.size(); ++i) {
			vzlib_out[i].init(true);
		}
	}

	// prepare Readstates OUT
	vstate_out.resize(num_split_files);

	unsigned seqlen = 0;
	std::string readstr;

	// loop until EOF - get reads from input files - write into split files
	for (auto inext = 0, iout = 0; next(inext, readstr, seqlen, true, orig_files);)
	{
		++vstate_out[iout].read_count;
		// if original files zipped, zip the splits too
		if (orig_files[inext].isZip) {
			int ret = vzlib_out[iout].defstr(readstr, ofsv[iout], vstate_out[iout].read_count == split_files[iout].numreads); // Z_STREAM_END | Z_OK - ok
			if (ret < Z_OK || ret > Z_STREAM_END) {
				ERR("Failed deflating readstring: ", readstr, " Output file idx: ", iout, " zlib status: ", ret);
				retval = false;
				break;
			}
		}
		else {
			ofsv[iout] << readstr;
			if (readstr.back() != '\n')
				ofsv[iout] << '\n';
		}

		// switch files
		if (is_two_files) inext ^= 1;
		if (vstate_out[iout].read_count == split_files[iout].numreads 
			&& ((is_paired && (iout & 1) == 1) || !is_paired)) 
			iout += num_sense; // switch split if previous split is done
		if (is_paired) iout ^= 1;
	} // ~for

	// close IN file streams
	for (unsigned i = 0; i < num_orig_files; ++i) {
		if (ifsv[i].is_open()) {
			ifsv[i].close();
		}
	}

	// close Out streams
	for (unsigned i = 0; i < ofsv.size(); ++i) {
		if (ofsv[i].is_open())
			ofsv[i].close();
	}

	// reset Readstates. No need for izlib - already cleaned
	for (unsigned i = 0; i < num_orig_files; ++i) {
		vstate_in[i].reset();
	}

	// zlib deflate streams already cleaned by now

	// calculate split file sizes
	for (unsigned i = 0; i < split_files.size(); ++i) {
		split_files[i].size = filesize(split_files[i].path.generic_string());
	}

	// write readfead descriptor
	write_descriptor();

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("Done splitting. Reads count: ", num_reads_tot, " Runtime sec: ", elapsed.count(), "\n");
	return retval;
} // ~Readfeed::split

// ---------------------------------------------------------------------------
// build_chunk_offsets  (INDEXED zipped)
//
// For each original gz file, do a single-pass decompression scan to collect
// the decompressed byte offset after every newline.  Then divide those line
// offsets into num_splits record-aligned chunks and store the byte boundaries
// in gz_slots[].
// ---------------------------------------------------------------------------
void Readfeed::build_chunk_offsets()
{
	auto start = std::chrono::high_resolution_clock::now();
	INFO("build_chunk_offsets: computing byte-range chunks for ", 
            num_orig_files, " file(s) x ", num_splits, " split(s)");

	const int linesPerRecord = orig_files[0].isFastq ? 4 : 2;
	// For a single interleaved paired file, chunk boundaries must align to complete pairs
	// (FWD + REV), so FWD and REV slot readers stay in sync when they share a reader.
	const bool is_interleaved = (num_orig_files < num_sense);
	const int alignUnit = is_interleaved ? linesPerRecord * static_cast<int>(num_sense) : linesPerRecord;

	gz_slots.resize(num_split_files);
	gz_slot_files.resize(num_split_files);

	for (size_t j = 0; j < num_orig_files; ++j) {
		auto& origFile = orig_files[j];

		// Pre-populate slot metadata for this file's sense (FWD slot only for interleaved)
		for (size_t i = 0; i < num_splits; ++i) {
			size_t slotIdx = i * num_sense + j;
			gz_slot_files[slotIdx].path     = origFile.path;
			gz_slot_files[slotIdx].isZip    = true;
			gz_slot_files[slotIdx].isFastq  = origFile.isFastq;
			gz_slot_files[slotIdx].isFasta  = origFile.isFasta;
			gz_slots[slotIdx].file_path     = origFile.path.generic_string();
		}

		// Pass 1: decompress entire file, collect newline offsets
		INFO("scanning ", origFile.path.generic_string());
		std::vector<uint64_t> newlineEnds;
		newlineEnds.reserve(static_cast<size_t>(origFile.numreads) * linesPerRecord + 1);

		{
			using PGR = rapidgzip::ParallelGzipReader<>;
			auto reader = std::make_unique<PGR>(
				std::make_unique<rapidgzip::StandardFileReader>(origFile.path.generic_string()),
				static_cast<size_t>(num_splits)
			);
			constexpr size_t CHUNK = 1U << 20; // 1 MiB
			std::vector<uint8_t> buf(CHUNK);
			uint64_t pos = 0;
			for (;;) {
				auto n = reader->read(reinterpret_cast<char*>(buf.data()), CHUNK);
				if (n <= 0) break;
				for (size_t k = 0; k < static_cast<size_t>(n); ++k) {
					if (buf[k] == '\n') newlineEnds.push_back(pos + k + 1);
				}
				pos += static_cast<uint64_t>(n);
			}
		}

		const uint64_t totalLines = static_cast<uint64_t>(newlineEnds.size());
		INFO("build_chunk_offsets: file ", j, " totalLines=", totalLines);

		// Compute line boundaries for each split, aligned to alignUnit (pair for interleaved)
		std::vector<uint64_t> alignedStarts(num_splits + 1);
		alignedStarts[0] = 0;
		for (size_t i = 1; i < num_splits; ++i) {
			uint64_t nominalLine = (totalLines * i) / num_splits;
			alignedStarts[i] = (nominalLine / alignUnit) * alignUnit;
		}
		alignedStarts[num_splits] = totalLines;

		for (size_t i = 0; i < num_splits; ++i) {
			size_t slotIdx = i * num_sense + j;
			uint64_t startLine = alignedStarts[i];
			uint64_t endLine   = alignedStarts[i + 1];

			uint64_t startByte = (startLine == 0 || newlineEnds.empty()) ? 0 : newlineEnds[startLine - 1];
			uint64_t endByte   = (endLine   == 0 || newlineEnds.empty()) ? 0 : newlineEnds[endLine   - 1];

			gz_slots[slotIdx].bytes_start = startByte;
			gz_slots[slotIdx].bytes_end   = endByte;
			gz_slot_files[slotIdx].numreads = static_cast<unsigned>((endLine - startLine) / linesPerRecord);

			INFO("build_chunk_offsets: slot ", slotIdx,
				" bytes=[", startByte, ",", endByte, ")",
				" reads=", gz_slot_files[slotIdx].numreads);
		}
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("build_chunk_offsets done in ", elapsed.count(), " sec");
} // ~Readfeed::build_chunk_offsets

// ---------------------------------------------------------------------------
// build_flat_chunk_offsets  (INDEXED_FLAT)
//
// For each original flat file, do a single-pass byte scan to collect the byte
// offset after every newline.  Then divide those line offsets into num_splits
// record-aligned chunks and store the byte boundaries in flat_slots[].
// ---------------------------------------------------------------------------
void Readfeed::build_flat_chunk_offsets()
{
	auto start = std::chrono::high_resolution_clock::now();
	INFO("build_flat_chunk_offsets: computing byte-range chunks for ", num_orig_files, 
            " file(s) x ", num_splits, " split(s)", " num_split_files | number of slots = ", num_split_files);

	const int linesPerRecord = orig_files[0].isFastq ? 4 : 2;
	// For a single interleaved paired file, chunk boundaries must align to complete pairs
	// (FWD + REV), so FWD and REV slot readers stay in sync when they share a reader.
	const bool is_interleaved = (num_orig_files < num_sense);
	const int alignUnit = is_interleaved ? linesPerRecord * static_cast<int>(num_sense) : linesPerRecord;

	flat_slots.resize(num_split_files);
	flat_slot_files.resize(num_split_files);

	for (size_t j = 0; j < num_orig_files; ++j) {
		auto& origFile = orig_files[j];

		// Pre-populate slot metadata for this file's sense (FWD slot only for interleaved)
		for (size_t i = 0; i < num_splits; ++i) {
			size_t slotIdx = i * num_sense + j;
			flat_slot_files[slotIdx].path    = origFile.path;
			flat_slot_files[slotIdx].isZip   = false;
			flat_slot_files[slotIdx].isFastq = origFile.isFastq;
			flat_slot_files[slotIdx].isFasta = origFile.isFasta;
			flat_slots[slotIdx].file_path    = origFile.path.generic_string();
		}

		// Pass 1: scan file, collect newline offsets
		INFO("build_flat_chunk_offsets: scanning ", origFile.path.generic_string());
		std::vector<uint64_t> newlineEnds;
		newlineEnds.reserve(static_cast<size_t>(origFile.numreads) * linesPerRecord + 1);

		{
			std::ifstream ifs(origFile.path, std::ios_base::in | std::ios_base::binary);
			if (!ifs.is_open()) {
				ERR("failed to open: ", origFile.path.generic_string());
				exit(1);
			}
			constexpr size_t CHUNK = 1U << 20; // 1 MiB
			std::vector<char> buf(CHUNK);
			uint64_t pos = 0;
			for (;;) {
				ifs.read(buf.data(), static_cast<std::streamsize>(CHUNK));
				auto n = ifs.gcount();
				if (n <= 0) break;
				for (size_t k = 0; k < static_cast<size_t>(n); ++k) {
					if (buf[k] == '\n') newlineEnds.push_back(pos + k + 1);
				}
				pos += static_cast<uint64_t>(n);
			}
		}

		const uint64_t totalLines = static_cast<uint64_t>(newlineEnds.size());
		INFO("build_flat_chunk_offsets: file ", j, " totalLines=", totalLines);

		// Compute line boundaries for each split, aligned to alignUnit (pair for interleaved)
		std::vector<uint64_t> alignedStarts(num_splits + 1);
		alignedStarts[0] = 0;
		for (size_t i = 1; i < num_splits; ++i) {
			uint64_t nominalLine = (totalLines * i) / num_splits;
			alignedStarts[i] = (nominalLine / alignUnit) * alignUnit;
		}
		alignedStarts[num_splits] = totalLines;

		for (size_t i = 0; i < num_splits; ++i) {
			size_t slotIdx = i * num_sense + j;
			uint64_t startLine = alignedStarts[i];
			uint64_t endLine   = alignedStarts[i + 1];

			uint64_t startByte = (startLine == 0 || newlineEnds.empty()) ? 0 : newlineEnds[startLine - 1];
			uint64_t endByte   = (endLine   == 0 || newlineEnds.empty()) ? 0 : newlineEnds[endLine   - 1];

			flat_slots[slotIdx].bytes_start = startByte;
			flat_slots[slotIdx].bytes_end   = endByte;
			flat_slot_files[slotIdx].numreads = static_cast<unsigned>((endLine - startLine) / linesPerRecord);

			INFO("build_flat_chunk_offsets: slot ", slotIdx,
				" bytes=[", startByte, ",", endByte, ")",
				" reads=", flat_slot_files[slotIdx].numreads);
		}
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("build_flat_chunk_offsets done in ", elapsed.count(), " sec");
} // ~Readfeed::build_flat_chunk_offsets

bool Readfeed::is_split_ready() {
	// compare with data in descriptor
	std::ifstream ifs; // descriptor file
	std::filesystem::path fn = basedir / "readfeed";
	if (std::filesystem::exists(fn)) {
		INFO("found existing readfeed descriptor ", fn.generic_string());
		ifs.open(fn, std::ios_base::in | std::ios_base::binary);
		if (!ifs.is_open()) {
			ERR("failed to open: ", fn.generic_string());
			exit(1);
		}

		unsigned lidx = 0; // line index
		unsigned fcnt = 0; // count of file entries in the descriptor
		unsigned fcnt_max = num_splits * num_sense + num_orig_files;
		unsigned fidx = 0; // file index
		unsigned fpidx = 0; // index of file parameters: name, size, lines, zip, fastq/fasta
		for (std::string line; std::getline(ifs, line); ) {
			if (!line.empty()) {
				// trim
				line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
				if (line[0] != '#') { // skip comments
					if (lidx == 0)
						is_ready = true; // timestamp
					else if (lidx == 1)
						is_ready = is_ready && std::stoul(line) == num_orig_files;
					else if (lidx == 2)
						is_ready = is_ready && std::stoul(line) == num_sense;
					else if (lidx == 3)
						is_ready = is_ready && std::stoul(line) == num_splits;
					else if (lidx == 4)
						is_ready = is_ready && std::stoul(line) == num_reads_tot;
					else {
						if (fpidx == 0) { // file path
							++fcnt;
							if (fcnt == num_orig_files+1) fidx = 0; // switch file index orig -> split files
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
								is_ready = is_ready && std::stoull(line) == orig_files[fidx].size; // original files
							else
								split_files[fidx].size = std::stoull(line); // set sizes of split files
						}
						else if (fpidx == 2) { // number of reads
							if (fcnt <= num_orig_files)
								is_ready = is_ready && std::stoul(line) == orig_files[fidx].numreads;
							else
								is_ready = is_ready && std::stoul(line) == split_files[fidx].numreads;
						}
						else if (fpidx == 3) { // is zip
							bool isZip = std::stoi(line) == 1;
							if (fcnt <= num_orig_files)
								is_ready = is_ready && orig_files[fidx].isZip == isZip;
							else
								is_ready = is_ready && split_files[fidx].isZip == isZip;
						}
						else if (fpidx == 4) { // fastq/a
							bool isFastq = line == "fastq";
							bool isFasta = line == "fasta";
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
			for (std::size_t i = 0; i < split_files.size(); ++i) {
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

bool Readfeed::define_format(const int& dbg)
{
	bool is_ascii = true;
	is_format_defined = true;

	ifsv.resize(num_orig_files);
	vstate_in.resize(num_orig_files);

	for (decltype(num_orig_files) i = 0; i < num_orig_files; ++i) {
		if (!ifsv[i].is_open()) {
			ifsv[i].open(orig_files[i].path, std::ios_base::in | std::ios_base::binary);
		}
		if (!ifsv[i].is_open()) {
			ERR("Failed to open file ", orig_files[i].path);
			exit(1);
		}
		auto fsz = std::filesystem::file_size(orig_files[i].path);
		auto blen = fsz > 100 ? 100 : fsz; // num bytes to read: max 100 - issue 290  20210511
		std::string str(100, '\0');
		ifsv[i].read(&str[0], blen); // get blen bytes from the stream
		if (dbg > 1) {
			auto st = ifsv[i].rdstate();
			INFO("rdstate: ", st); // 3 - some undefined state. Defined ones are 0,1,2,4
		}
		for (std::size_t i = 0; i < str.size(); ++i) {
			// 20201008 TODO: this is quite adhoc - need a better validation like evaluating the gz, zlib header
			// warning: comparison is always false due to limited range of data type [-Wtype-limits]
			if (str[i] > 127 || str[i] < 0) {
				is_ascii = false; // if any of the char is not ascii -> file is not ascii
				break;
			}
		}
		orig_files[i].isZip = !is_ascii;
		if (!is_ascii) {
			// init izlib for inflation
			if (vzlib_in.size() < (size_t)i+1) {
				vzlib_in.emplace_back(Izlib());
			}
			vzlib_in[i].init(false);
		}
		
		ifsv[i].seekg(0); // rewind stream to the start
		// get a read to test ability to read
		if (next(i, str, true, orig_files)) {
			std::string fmt = orig_files[i].isZip ? "gzipped" : "flat ASCII";
			if (orig_files[i].isFasta) {
				//biof = BIO_FORMAT::FASTA;
				INFO("file: ", orig_files[i].path, " is FASTA ", fmt);
			}
			else if (orig_files[i].isFastq) {
				//biof = BIO_FORMAT::FASTQ;
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
			vzlib_in[i].reset_inflate();
		}
	} // ~for orig files

	return is_format_defined;
} // ~Readfeed::define_format

// ---------------------------------------------------------------------------
// count_reads_parallel  (FEED_TYPE::INDEXED)
//
// For gzipped input files: creates a ParallelGzipReader per file with
// num_splits decompression workers and scans the decompressed stream in a
// single pass to count records and accumulate sequence lengths.
//
// For flat input files: finds record-aligned byte boundaries for each split,
// then spawns num_splits threads that each scan their assigned byte range.
//
// Populates: num_reads_tot, length_all, min_read_len, max_read_len, and
//            orig_files[j].numreads.
// ---------------------------------------------------------------------------
void Readfeed::count_reads_parallel()
{
	auto start = std::chrono::high_resolution_clock::now();
	INFO("count_reads_parallel: started with ", num_splits, " worker thread(s)");

	if (!is_format_defined)
		define_format();

	const int linesPerRecord = orig_files[0].isFastq ? 4 : 2;

	for (size_t j = 0; j < num_orig_files; ++j) {
		auto& origFile = orig_files[j];

		if (origFile.isZip) {
			// --- gz: ParallelGzipReader decompresses with num_splits threads internally ---
			using PGR = rapidgzip::ParallelGzipReader<>;
			auto reader = std::make_unique<PGR>(
				std::make_unique<rapidgzip::StandardFileReader>(origFile.path.generic_string()),
				static_cast<size_t>(num_splits)
			);

			constexpr size_t CHUNK = 1U << 20; // 1 MiB
			std::vector<uint8_t> buf(CHUNK);

			int lineInRecord = 0; // 0=header, 1=seq, 2='+', 3=qual (FASTQ); 0=header, 1=seq (FASTA)
			uint64_t seqlen = 0;

			for (;;) {
				const auto n = reader->read(reinterpret_cast<char*>(buf.data()), CHUNK);
				if (n <= 0) break;
				for (size_t k = 0; k < static_cast<size_t>(n); ++k) {
					if (buf[k] == '\n') {
						if (lineInRecord == 1) { // end of sequence line
							++origFile.numreads;
							length_all += seqlen;
							if (seqlen > max_read_len) max_read_len = static_cast<uint32_t>(seqlen);
							if (min_read_len == 0 || static_cast<uint32_t>(seqlen) < min_read_len)
								min_read_len = static_cast<uint32_t>(seqlen);
							seqlen = 0;
						}
						lineInRecord = (lineInRecord + 1) % linesPerRecord;
					} else if (lineInRecord == 1) {
						++seqlen;
					}
				}
			}
		} else {
			// --- flat: parallel threads, each scanning a record-aligned byte range ---
			const uint64_t fileSize = static_cast<uint64_t>(origFile.size);

			// Build record-aligned split boundaries.
			// boundaries[i]   = byte offset where thread i starts (inclusive)
			// boundaries[i+1] = byte offset where thread i ends   (exclusive)
			std::vector<uint64_t> boundaries(num_splits + 1);
			boundaries[0] = 0;
			boundaries[num_splits] = fileSize;

			if (num_splits > 1) {
				// Sequential boundary-finding pass: seek near each split point and advance
				// to the start of the next complete record.
				std::ifstream bifs(origFile.path, std::ios_base::in | std::ios_base::binary);
				if (!bifs.is_open()) {
					ERR("count_reads_parallel: cannot open ", origFile.path.generic_string());
					exit(1);
				}

				for (size_t i = 1; i < num_splits; ++i) {
					bifs.seekg(static_cast<std::streamoff>(fileSize * i / num_splits));
					std::string ln;
					std::getline(bifs, ln); // skip to end of current (partial) line

					if (origFile.isFastq) {
						// Read 4-line groups until we find one where line[0] starts with '@'
						// and line[2] starts with '+' (valid FASTQ record start).
						bool found = false;
						for (int attempt = 0; attempt < linesPerRecord && !found; ++attempt) {
							const uint64_t candidatePos = static_cast<uint64_t>(bifs.tellg());
							std::string l0, l1, l2, l3;
							if (!std::getline(bifs, l0)) { boundaries[i] = fileSize; break; }
							if (!std::getline(bifs, l1)) { boundaries[i] = fileSize; break; }
							if (!std::getline(bifs, l2)) { boundaries[i] = fileSize; break; }
							if (!std::getline(bifs, l3)) { boundaries[i] = fileSize; break; }
							if (!l0.empty() && l0[0] == '@' && !l2.empty() && l2[0] == '+') {
								boundaries[i] = candidatePos;
								found = true;
							} else {
								// Advance by one line and retry
								bifs.seekg(static_cast<std::streamoff>(candidatePos + l0.size() + 1));
							}
						}
						if (!found) boundaries[i] = fileSize; // fold empty range into previous split
					} else {
						// FASTA: scan for the next '>' at the start of a line
						bool found = false;
						std::string fln;
						while (!found) {
							const uint64_t pos = static_cast<uint64_t>(bifs.tellg());
							if (!std::getline(bifs, fln)) { boundaries[i] = fileSize; break; }
							if (!fln.empty() && fln[0] == '>') { boundaries[i] = pos; found = true; }
						}
					}
				}
			}

			// Thread-local accumulator
			struct ThreadResult {
				uint64_t numreads = 0;
				uint64_t length   = 0;
				uint32_t minlen   = 0;
				uint32_t maxlen   = 0;
			};
			std::vector<ThreadResult> results(num_splits);

			{
				std::vector<std::thread> workers;
				workers.reserve(num_splits);
				for (size_t i = 0; i < num_splits; ++i) {
					workers.emplace_back([&, i]() {
						auto& res = results[i];
						const uint64_t startByte = boundaries[i];
						const uint64_t endByte   = boundaries[i + 1];
						if (startByte >= endByte) return;

						std::ifstream ifs(origFile.path, std::ios_base::in | std::ios_base::binary);
						if (!ifs.is_open()) return;
						ifs.seekg(static_cast<std::streamoff>(startByte));

						constexpr size_t BUFSZ = 1U << 16; // 64 KiB
						std::vector<char> buf(BUFSZ);
						uint64_t remaining = endByte - startByte;
						int lineInRecord = 0;
						uint64_t seqlen = 0;

						while (remaining > 0) {
							const size_t toRead = static_cast<size_t>(
								std::min(static_cast<uint64_t>(BUFSZ), remaining));
							ifs.read(buf.data(), static_cast<std::streamsize>(toRead));
							const auto n = static_cast<size_t>(ifs.gcount());
							if (n == 0) break;
							remaining -= n;
							for (size_t k = 0; k < n; ++k) {
								if (buf[k] == '\n') {
									if (lineInRecord == 1) { // end of sequence line
										++res.numreads;
										res.length += seqlen;
										if (seqlen > res.maxlen) res.maxlen = static_cast<uint32_t>(seqlen);
										if (res.minlen == 0 || static_cast<uint32_t>(seqlen) < res.minlen)
											res.minlen = static_cast<uint32_t>(seqlen);
										seqlen = 0;
									}
									lineInRecord = (lineInRecord + 1) % linesPerRecord;
								} else if (lineInRecord == 1) {
									++seqlen;
								}
							}
						}
					});
				}
				for (auto& t : workers) t.join();
			}

			// Reduce per-thread results into origFile and class-level members
			for (const auto& r : results) {
				origFile.numreads += static_cast<unsigned>(r.numreads);
				length_all += r.length;
				if (r.maxlen > max_read_len) max_read_len = r.maxlen;
				if (min_read_len == 0 || (r.minlen > 0 && r.minlen < min_read_len))
					min_read_len = r.minlen;
			}
		}

		num_reads_tot += origFile.numreads;
	}

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("count_reads_parallel done. Elapsed: ", elapsed.count(),
	     " sec. Total reads: ", num_reads_tot);
} // ~Readfeed::count_reads_parallel

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
	for (std::size_t i = 0; i < orig_files.size(); ++i) {
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
		if (is_two_files) inext ^= 1; // toggle the index of the input file
	} // ~for

	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("done count. Elapsed time: ", elapsed.count(), " sec. Total reads: ", num_reads_tot);

} // ~Readfeed::count_reads

/*
  init split files' names, and attributes including format fastq/a, zip, number of reads.
*/
void Readfeed::init_split_files()
{
	split_files.resize(num_split_files);
	std::stringstream ss;

	for (decltype(num_splits) i = 0, idx = 0; i < num_splits; ++i) {
		for (decltype(num_sense) j = 0; j < num_sense; ++j, ++idx) {
			std::string stem = j == 0 ? "fwd_" : "rev_";
			auto jj = is_two_files ? j : 0;
			std::string sfx_1 = orig_files[jj].isFasta ? ".fa" : ".fq";
			auto sfx = orig_files[jj].isZip ? sfx_1 + ".gz" : sfx_1;
			ss << stem << i << sfx; // split file basename
			//auto idx = i * num_orig_files + j;
			split_files[idx].path = basedir / ss.str();
			ss.str("");
			INFO("adding split to the list: ", split_files[idx].path.generic_string());
		}
	}

	// calculate number of reads in each of the split files
	auto nreads = num_reads_tot / num_sense; // num reads of the same sense e.g. FWD
	auto minr = nreads / num_splits; // quotient i.e. min number of reads in each output file
	auto surplus = nreads - minr * num_splits; // remainder of reads to be distributed between the output files
    // for each split
	for (size_t i = 0, idx = 0; i < num_splits; ++i) {
		auto maxr = i < surplus ? minr + 1 : minr; // distribute the surplus
        // for each sense
		for (size_t j = 0; j < num_sense; ++j, ++idx) {
			split_files[idx].numreads = maxr;
			auto jj = is_two_files ? j : 0;
			split_files[idx].isZip = orig_files[jj].isZip;
			split_files[idx].isFastq = orig_files[jj].isFastq;
			split_files[idx].isFasta = orig_files[jj].isFasta;
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
		"#   num_sense\n"
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
		std::ofstream ofs(fn, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!ofs.is_open()) {
			ERR("failed to open file: ", fn.generic_string());
			exit(1);
		}

        INFO("writing reads descriptor to: ", fn.generic_string());
		std::time_t tm = std::time(0);
		ofs << comments << '\n';
		ofs << std::ctime(&tm) << '\n'; // line 0: timestamp  <ctime> 'Tue Oct 20 08:39:35 2020'
		ofs << num_orig_files << '\n';  // line 1:
		ofs << num_sense << '\n';       // line 2
		ofs << num_splits << '\n';      // line 3:
		ofs << num_reads_tot << '\n';   // line 4:
		for (std::size_t i = 0; i < orig_files.size(); ++i) {
			ofs << orig_files[i].path.generic_string() << '\n';     // file path
			ofs << orig_files[i].size << '\n';     // file size
			ofs << orig_files[i].numreads << '\n'; // reads count
			ofs << orig_files[i].isZip << '\n';    // zip
			auto fx = orig_files[i].isFastq == true ? "fastq" : "fasta";
			ofs << fx << '\n';
		}
		for (std::size_t i = 0; i < split_files.size(); ++i) {
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
	if (type == FEED_TYPE::INDEXED) return;

	vzlib_in.resize(split_files.size());
	for (std::size_t i = 0; i < vzlib_in.size(); ++i) {
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
	if (type == FEED_TYPE::INDEXED && orig_files[0].isZip) {
		vstate_in.resize(gz_slots.size());
		for (auto& s : vstate_in) s.reset();

		// For interleaved paired, REV slots (odd) share the FWD slot's reader — skip them.
		const bool is_interleaved = (num_orig_files < num_sense);
		for (std::size_t i = 0; i < gz_slots.size(); ++i) {
			if (is_interleaved && i % num_sense != 0) continue;
			auto& slot = gz_slots[i];
			slot.reader = std::unique_ptr<GzReaderImpl, GzReaderDeleter>(
				new GzReaderImpl(
					std::make_unique<rapidgzip::StandardFileReader>(slot.file_path),
					/*parallelization=*/std::size_t(1)
				)
			);
			slot.reader->rdr.seek(static_cast<long long>(slot.bytes_start));
			slot.bytes_remaining = slot.bytes_end - slot.bytes_start;
			slot.buf_pos = 0;
			slot.buf_len = 0;
		}
		return;
	}

	if (type == FEED_TYPE::INDEXED && orig_files[0].isZip == false) {
		vstate_in.resize(flat_slots.size());
		for (auto& s : vstate_in) s.reset();

		// For interleaved paired, REV slots (odd) share the FWD slot's ifstream — skip them.
		const bool is_interleaved = (num_orig_files < num_sense);
		for (std::size_t i = 0; i < flat_slots.size(); ++i) {
			if (is_interleaved && i % num_sense != 0) continue;
			auto& slot = flat_slots[i];
			if (slot.ifs.is_open()) slot.ifs.close();
			slot.ifs.open(slot.file_path, std::ios_base::in | std::ios_base::binary);
			if (!slot.ifs.is_open()) {
				ERR("failed to open: ", slot.file_path);
				exit(1);
			}
			slot.ifs.seekg(static_cast<std::streamoff>(slot.bytes_start));
			slot.bytes_remaining = slot.bytes_end - slot.bytes_start;
			slot.buf_pos = 0;
			slot.buf_len = 0;
		}
		return;
	}

    // split files - outdated - to be removed
	for (std::size_t i = 0; i < vstate_in.size(); ++i) {
		vstate_in[i].reset();
	}
	vstate_in.resize(split_files.size());

	// ifsv
	for (std::size_t i = 0; i < ifsv.size(); ++i) {
		if (ifsv[i].is_open()) ifsv[i].close();
	}
	ifsv.resize(split_files.size());
	for (std::size_t i = 0; i < split_files.size(); ++i) {
		ifsv[i].open(split_files[i].path, std::ios_base::in | std::ios_base::binary);
		if (!ifsv[i].is_open()) {
			ERR("failed to open: ", split_files[i].path.generic_string());
			exit(1);
		}
	}

	// vzlib_in
	init_vzlib_in();
} // ~Readfeed::init_reading

int Readfeed::clean()
{
	std::ifstream ifs; // descriptor file
	std::filesystem::path fn = basedir / "readfeed";
	int n_del = 0; // count of deleted files
	if (std::filesystem::exists(fn)) {
		INFO("found descriptor ", fn.generic_string());
		ifs.open(fn, std::ios_base::in | std::ios_base::binary);
		if (!ifs.is_open()) {
			ERR("failed to open: ", fn.generic_string());
			exit(1);
		}

		int lidx = 0; // line index
		int fcnt = 0; // count of file entries in the desctiptor
		int fcnt_max = 0; // number of All file entries in the descriptor
		int fidx = 0; // file index
		int fpidx = 0; // index of file parameters: name, size, lines, zip, fastq/fasta
		int n_orig = 0; // number of original files
		int n_sense = 0; // number of senses
		int n_split = 0; // number of splits
		//unsigned n_tot = 0; // total of reads in orig files
		std::regex rx_num("[0-9]+"); // (-|+)|][0-9]+

		for (std::string line; std::getline(ifs, line); ) {
			if (!line.empty()) {
				// trim
				line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
				if (line[0] != '#') { // skip comments
					if (lidx == 1 || lidx == 2 || lidx == 3 || lidx == 4) {
						if (!std::regex_match(line, rx_num)) {
							ERR("not a number: '", line, "'");
							exit(1);
						}
					}
					if (lidx == 0); // skip timestamp
					else if (lidx == 1) n_orig = std::stoi(line);
					else if (lidx == 2) n_sense = std::stoi(line);
					else if (lidx == 3) n_split = std::stoi(line);
					else if (lidx == 4) {
						//n_tot = std::stoi(line);
						fcnt_max = n_split * n_sense + n_orig;
					}
					else {
						if (fpidx == 0) { // file path
							++fcnt;
							if (fcnt > n_orig) {
								// delete a split file
								std::filesystem::path sf(line);
								if (std::filesystem::exists(sf) && sf.parent_path() == basedir) {
									INFO("removing split file: ", line);
									std::filesystem::remove(sf);
									++n_del;
								}
							}
							if (fcnt == n_orig + 1)	fidx = 0; // switch file index orig -> split files
							if (fcnt > fcnt_max) {
								INFO("expected max file count: ", fcnt_max, " current file count in the descriptor: ", fcnt);
								is_ready = false;
								break;
							}
						}
						else if (fpidx == 4) { // fastq/a
							fpidx = 0;
							++fidx;
							++lidx;
							continue;
						}
						++fpidx;
					}
					++lidx;
				}
			}
		} // ~for lines in descriptor
	}
	return n_del;
} // ~Readfeed::clean
