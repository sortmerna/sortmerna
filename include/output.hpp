#pragma once
/**
 * FILE: output.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Handles output to files
 */
#include <fstream>
#include <stdint.h>
#include <string>

#include "common.hpp"
#include "options.hpp"
#include "readstats.hpp"
//#include "read.hpp"

// forward
struct Index;
class References;
class Read;

class Output {
public:
	// output streams for aligned reads (FASTA/FASTQ, SAM and BLAST-like)
	ofstream acceptedreads;
	ofstream acceptedsam;
	ofstream acceptedblast;
	ofstream logout;

	// file names
	std::string acceptedstrings;
	std::string acceptedstrings_sam; // used in Index::load_stats
	std::string acceptedstrings_blast;
	std::string logoutfile;
	std::string denovo_otus_file;
	std::string acceptedotumap_file;

	Output(Runopts & opts, Readstats & readstats)
		:
		opts(opts)
	{
		init(readstats);
	}
	~Output() {}

	void init(Readstats & readstats); // TODO: make private?

	void report_blast(
		ofstream & fileout,
		Index & index,
		References & refs,
		Read & read
	);

	void report_sam(
		ofstream & fileout,
		References & refs,
		Read & read
	);

	void calcMismatchGapId(References & refs, Read & read, int alignIdx, uint32_t & mismatches, uint32_t & gaps, double & id);
	void openfiles();
	void closefiles();

private:
	Runopts & opts;
}; // ~class Output


void generateReports(Runopts opts);