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
		bool isFastq = true;
		bool lastRec = false; // lastRec is to make one iteration past the EOF
		Gzip gzip(opts);

		ss << id << " thread: " << std::this_thread::get_id() << " started\n";
		std::cout << ss.str(); ss.str("");
		auto t = std::chrono::high_resolution_clock::now();

		// read lines from the files and create read objects
		for (int count = 0, stat = 0; ; ) // count lines in a single record
		{
			stat = gzip.getline(ifs, line);
			if (stat == RL_END)
			{
				// push the last ready read object to the queue
				if (!read.isEmpty)
				{
					read.init(opts, kvdb, read_id); // load alignment statistics from DB
					readQueue.push(read);
					++read_id;
				}
				break;
			}

			if (stat == RL_ERR)
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " ERROR reading from Reads file. Exiting..." << std::endl;
				exit(1);
			}

			if (line.empty()) continue;

			// left trim space and '>' or '@'
			//line.erase(line.begin(), std::find_if(line.begin(), line.end(), [](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));
			// right-trim whitespace in place (removes '\r' too)
			line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
			// removes all space
			//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));

			// fastq: 0(header), 1(seq), 2(+), 3(quality)
			// fasta: 0(header), 1(seq)
			if (line[0] == FASTA_HEADER_START || line[0] == FASTQ_HEADER_START)
			{ // add header -->
				if (!read.isEmpty)
				{ // push previous read object to queue
					read.init(opts, kvdb, read_id);
					readQueue.push(read);
					++read_id;
					count = 0;
				}

				// start new record
				read.clear();
				isFastq = (line[0] == FASTQ_HEADER_START);
				read.format = isFastq ? Format::FASTQ : Format::FASTA;
				read.header = line;
				read.isEmpty = false;
			} // ~if header line
			else 
			{ // add sequence -->
				++count;
				if (isFastq && line[0] == '+') continue;
				if (isFastq && count == 3)
				{
					read.quality = line;
					continue;
				}
				read.sequence += line;
			}
			//if (ifs.eof()) lastRec = true; // push and break
		} // ~for getline

		std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
		readQueue.mDoneAdding();
		ss << id << " thread: " << std::this_thread::get_id() << " done. Elapsed time: " 
			<< std::setprecision(2) << std::fixed << elapsed.count() << " sec Reads added: " << read_id << std::endl;
		std::cout << ss.str(); ss.str("");
	}
	ifs.close();
} // ~Reader::read