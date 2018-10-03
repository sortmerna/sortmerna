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

#include "reader.hpp"
#include "gzip.hpp"

void Reader::read()
{
	std::stringstream ss;

	std::ifstream ifs(opts.readsfile, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) {
		std::cerr << __FILE__ << ":" << __LINE__ << " failed to open " << opts.readsfile << std::endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		std::string line;
		Read read; // an empty read
		unsigned int read_id = 0; // read ID
		unsigned int tcount = 0;
		bool isFastq = false;
		bool isFasta = false;
		//bool lastRec = false; // lastRec is to make one iteration past the EOF
		Gzip gzip(opts); // reads both zipped and non-zipped files

		ss << id << " thread: " << std::this_thread::get_id() << " started\n";
		std::cout << ss.str(); ss.str("");
		auto t = std::chrono::high_resolution_clock::now();

		// read lines from the files and create read objects
		// NOTE: don't increment count here to avoid counting (just in case) empty lines
		for (int count = 0, stat = 0; ;	++count) // count lines in a single record
		{
			stat = gzip.getline(ifs, line);
			++tcount;

			if (stat == RL_END)
			{
				// push the last Read to the queue
				if (!read.isEmpty)
				{
					read.init(opts, kvdb, read_id); // load alignment statistics from DB
					readQueue.push(read);
				}
				break;
			}

			if (stat == RL_ERR)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " ERROR reading from Reads file. Exiting..." << std::endl;
				exit(1);
			}

			if (line.empty()) 
			{
				--count;
				--tcount;
				continue;
			}

			// left trim space and '>' or '@'
			//line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));
			// right-trim whitespace in place (removes '\r' too)
			line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
			// removes all space
			//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));
			if (tcount == 1)
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
			if ((isFasta && line[0] == FASTA_HEADER_START) || (isFastq && count == 0))
			{ // add header -->
				if (!read.isEmpty)
				{ // push previous read object to queue
					read.init(opts, kvdb, read_id);
					readQueue.push(read);
					++read_id;
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
					if (count == 2) // line[0] == '+' validation is already by readstats::calculate
						continue;
					if (count == 3)
					{
						read.quality = line;
						continue;
					}
				}

				read.sequence += line; // FASTA multi-line sequence or FASTQ sequence
			}
			//if (ifs.eof()) lastRec = true; // push and break
		} // ~for getline

		std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
		--readQueue.numPushers; // signal the reader done adding
		ss << id << " thread: " << std::this_thread::get_id() << " done. Elapsed time: " 
			<< std::setprecision(2) << std::fixed << elapsed.count() << " sec Reads added: " << read_id + 1 
			<< " readQueue.size: " << readQueue.size() << std::endl;
		std::cout << ss.str(); ss.str("");
	}
	ifs.close();
} // ~Reader::read

bool Reader::loadReadByIdx(Runopts & opts, Read & read)
{
	std::stringstream ss;
	bool isok = false;

	std::ifstream ifs(opts.readsfile, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) 
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " failed to open " << opts.readsfile << std::endl;
		exit(EXIT_FAILURE);
	}
	else
	{
		std::string line;
		unsigned int read_id = 0; // read ID
		bool isFastq = true;
		Gzip gzip(opts);

		auto t = std::chrono::high_resolution_clock::now();

		// read lines from the reads file
		for (int count = 0, stat = 0; ; ) // count lines in a single read
		{
			stat = gzip.getline(ifs, line);
			if (stat == RL_END) break;

			if (stat == RL_ERR)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " ERROR reading from Reads file. Exiting..." << std::endl;
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
				if (read_id == read.id)
				{
					isFastq = (line[0] == FASTQ_HEADER_START);
					read.format = isFastq ? Format::FASTQ : Format::FASTA;
					read.header = line;
					read.isEmpty = false;
				}
				else {
					++read_id;
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