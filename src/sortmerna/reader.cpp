/**
 * FILE: reader.cpp
 * Created: Nov 26, 2017 Sun
 */

#include <locale> // std::isspace

#include "reader.hpp"



void Reader::read()
{
	std::ifstream ifs(readsfile, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) {
		printf("failed to open %s\n", readsfile.c_str());
		exit(EXIT_FAILURE);
	}
	else
	{
		std::string line;
		Read read; // an empty read
		int id = 0; // read ID
		bool isFastq = true;
		bool lastRec = false; // lastRec is to make one iteration past the EOF

		std::stringstream ss;
		ss << std::this_thread::get_id();
		printf("Reader thread %s started\n", ss.str().c_str());
		ss.str(""); // clear the stream
		auto t = std::chrono::high_resolution_clock::now();

		for (int count = 0; ; ) // count lines in a single record
		{
			if (!lastRec) std::getline(ifs, line);

			if (line.empty() && !lastRec)
			{
				if (ifs.eof()) lastRec = true;
				continue;
			}

			if (lastRec)
			{
				if (!read.isEmpty)
				{
					read.init(opts, kvdb); // load alignment statistics from DB
					readQueue.push(read);
					++id;
				}
				break;
			}

			// remove whitespace in place (removes '\r' too)
			line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));
			// fastq: 0(header), 1(seq), 2(+), 3(quality)
			// fasta: 0(header), 1(seq)
			if (line[0] == FASTA_HEADER_START || line[0] == FASTQ_HEADER_START)
			{
				if (!read.isEmpty)
				{
					read.init(opts, kvdb);
					readQueue.push(read);
					++id;
					count = 0;
				}

				// start new record
				read.clear();
				read.id = id;
				isFastq = (line[0] == FASTQ_HEADER_START);
				read.format = isFastq ? Format::FASTQ : Format::FASTA;
				read.header = line;
				read.isEmpty = false;
			} // ~if header line
			else {
				++count;
				if (isFastq && line[0] == '+') continue;
				if (isFastq && count == 3)
				{
					read.quality = line;
					continue;
				}
				read.sequence += line;
			}
			if (ifs.eof()) lastRec = true; // push and break
		} // ~for getline

		std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
		readQueue.mDoneAdding();
		ss << std::this_thread::get_id();
		printf("Reader thread %s done. Elapsed time: %.2f sec Reads added: %d\n",
			ss.str().c_str(), elapsed.count(), id + 1);
	}
	ifs.close();
} // ~Reader::read