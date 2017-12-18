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

	void report_blast(
		Index & index,
		References & refs,
		Read & read
	);

	void report_sam(
		References & refs,
		Read & read
	);

	void writeSamHeader();

	void report_fasta();
	void report_denovo();
	void report_biom();

	void openfiles();
	void closefiles();

private:
	void init(Readstats & readstats);

private:
	Runopts & opts;
}; // ~class Output


void generateReports(Runopts opts);