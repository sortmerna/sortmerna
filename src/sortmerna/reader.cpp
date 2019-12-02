/**
 * FILE: reader.cpp
 * Created: Nov 26, 2017 Sun
 * @copyright 2016-19 Clarity Genomics BVBA
 * 
 * Processes Reads file, creates Read objects, and pushes them to a queue for further pick−up by Processor
 */

#include <vector>
#include <sstream> // std::stringstream
#include <ios> // std::ios_base
#include <algorithm> // find, find_if
#include <chrono> // std::chrono
#include <iomanip> // std::precision
#include <locale> // std::isspace

#include "reader.hpp"
#include "gzip.hpp"

Reader::Reader(std::string id, bool is_gzipped)
	:
	id(id),
	is_gzipped(is_gzipped),
	gzip(is_gzipped),
	is_done(false),
	read_count(0),
	line_count(0),
	last_count(0),
	last_stat(0),
	isFastq(false),
	isFasta(false)
{} // ~Reader::Reader

Reader::~Reader() {}

bool Reader::loadReadByIdx(Runopts & opts, Read & read)
{
	std::stringstream ss;
	bool isok = false;

	std::ifstream ifs(opts.readfiles[read.readfile_idx], std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) 
	{
		std::cerr << STAMP << "failed to open " << opts.readfiles[read.readfile_idx] << std::endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		std::string line;
		std::size_t read_num = 0;
		bool isFastq = true;
		Gzip gzip(opts.is_gz);

		auto t = std::chrono::high_resolution_clock::now();

		// read lines from the reads file
		for (int count = 0, stat = 0; ; ) // count lines in a single read
		{
			stat = gzip.getline(ifs, line);
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
					isFastq = (line[0] == FASTQ_HEADER_START);
					read.format = isFastq ? Format::FASTQ : Format::FASTA;
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
				if ( isFastq )
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
	}

	ifs.close();

	return isok;

} // ~Reader::loadRead

bool Reader::loadReadById(Runopts & opts, Read & read)
{
	return true;
} // ~Reader::loadReadById


/** 
 * return next read from the reads file on each call 
 */
Read Reader::nextread(std::ifstream &ifs, const uint8_t readsfile_idx, Runopts & opts)
{
	std::string line;
	Read read; // an empty read

	// read lines from the reads file and create Read object
	for (int count = last_count, stat = last_stat; !is_done; ++count) // count lines in a single record/read
	{
		if (last_header.size() > 0)
		{
			// start new record
			read.clear();
			read.format = isFastq ? Format::FASTQ : Format::FASTA;
			read.header = last_header;
			read.isEmpty = false;
			last_header.clear();
		}

		stat = gzip.getline(ifs, line);

		if (stat == RL_END)
		{
			// push the last Read to the queue
			if (!read.isEmpty)
			{
				read.read_num = read_count;
				read.readfile_idx = readsfile_idx;
				read.generate_id();
				++read_count;
				is_done = true;
			}
			break;
		}

		if (stat == RL_ERR)
		{
			std::cerr << STAMP << "ERROR reading from file: [" << opts.readfiles[readsfile_idx] << "]. Exiting..." << std::endl;
			exit(1);
		}

		if (line.empty())
		{
			--count;
			continue;
		}

		++line_count;

		// left trim space and '>' or '@'
		//line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));
		// right-trim whitespace in place (removes '\r' too)
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
		// removes all space
		//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));
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
		if ((isFasta && line[0] == FASTA_HEADER_START) || (isFastq && count == 0)) // header line reached
		{
			if (!read.isEmpty)
			{
				read.read_num = read_count;
				read.readfile_idx = readsfile_idx;
				read.generate_id();
				++read_count;
				last_header = line;
				last_count = 1;
				last_stat = stat;
				break; // return the read here
			}

			// start new record
			read.clear();
			read.format = isFastq ? Format::FASTQ : Format::FASTA;
			read.header = line;
			read.isEmpty = false;

			count = 0; // FASTA record start
		} // ~if header line
		else
		{ // add sequence -->
			if (isFastq)
			{
				if (count == 2) // line[0] == '+' validation is already done by readstats::calculate
					continue;
				if (count == 3)
				{
					read.quality = line;
					continue;
				}
			}
			read.sequence += line; // FASTA multi-line sequence or FASTQ sequence
		}
	} // ~for getline
	return read;
} // ~Reader::nextread

/**
 * get a next read sequence from the reads file
 * @return true if record exists, else false
 */
bool Reader::nextread(std::ifstream& ifs, const std::string &readsfile, std::string &seq)
{
	bool has_seq = false;
	seq = ""; // ensure empty
	std::string line;

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
	read_count = 0;
	line_count = 0;
	last_count = 0;
	last_stat = 0;
	bool isFastq = false;
	bool isFasta = false;
	is_done = false;
}