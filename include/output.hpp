#pragma once
/**
 * FILE: output.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Handles output to files
 */
#include <fstream>

#include "options.hpp"
#include "readstats.hpp"

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

private:
	Runopts & opts;
}; // ~class Output
