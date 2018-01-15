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
#include <vector>

#include "common.hpp"

// forward
struct Index;
class References;
class Read;
class Refstats;
struct Readstats;
struct Runopts;

class Output {
public:
	// output streams for aligned reads (FASTA/FASTQ, SAM and BLAST-like)
	std::ofstream acceptedreads;
	std::ofstream acceptedsam;
	std::ofstream acceptedblast;
	std::ofstream logstream;
	std::ofstream denovoreads;
	std::ofstream otherreads;
	std::ofstream biomout;

	// file names
	std::string acceptedstrings;
	std::string acceptedstrings_sam; // used in Index::load_stats
	std::string acceptedstrings_blast;
	std::string logfile;
	std::string denovo_otus_file;
	std::string acceptedotumap_file;
	std::string biomfile;

	Output(Runopts & opts, Readstats & readstats)
	{
		init(opts, readstats);
	}
	~Output() {}

	void report_blast(
		Runopts & opts,
		Refstats & refstats,
		References & refs,
		Read & read
	);

	void report_sam(
		Runopts & opts,
		References & refs,
		Read & read
	);

	void writeSamHeader(Runopts & opts);

	void report_fasta(Runopts & opts, std::vector<Read> & reads);
	void report_denovo(Runopts & opts, std::vector<Read> & reads);
	void report_biom();

	void openfiles(Runopts & opts);
	void closefiles();

private:
	void init(Runopts & opts, Readstats & readstats);

}; // ~class Output


void generateReports(Runopts & opts);