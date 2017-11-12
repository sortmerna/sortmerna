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

	// file names
	std::string acceptedstrings;
	std::string acceptedstrings_sam; // used in Index::load_stats
	std::string acceptedstrings_blast;
	std::string logoutfile;
	std::string denovo_otus_file;
	std::string acceptedotumap_file;

	uint64_t total_reads_mapped = 0; // shared by Processor threads
	uint64_t total_reads_mapped_cov = 0; // shared   total number of reads mapped passing E-value threshold and %id and/or %query coverage thresholds
	std::vector<uint64_t> reads_matched_per_db; // total number of reads matched for each database

	Output(Runopts & opts, Readstats & readstats)
		:
		opts(opts),
		reads_matched_per_db(opts.indexfiles.size(), 0)
	{
		init(readstats);
	}
	~Output() {}

	void init(Readstats & readstats); // TODO: make private?

private:
	Runopts & opts;
}; // ~class Output
