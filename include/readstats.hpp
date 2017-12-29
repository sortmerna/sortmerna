#pragma once
/**
 * FILE: readstats.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Collective Statistics for all Reads. Encapsulates old 'compute_read_stats' logic and results
 */

#include <cstdint>
#include <string>
#include <vector>
#include <map>

#include "common.hpp"
#include "options.hpp"

struct Readstats {
	Runopts & opts;

	// TODO: synchronize - Computed in worker thread (per read) - shared
	uint32_t min_read_len; // length of the shortest Read in the Reads file. 'parallelTraversalJob'
	uint32_t max_read_len; // length of the longest Read in the Reads file. 'parallelTraversalJob'
	// total number of reads mapped passing E-value threshold. Computed in 'compute_lis_alignment' in a worker thread i.e. per read.
	uint64_t total_reads_mapped; // Set in 'compute_lis_alignment2'
	// total number of reads mapped passing E - value threshold and %id and/or %query coverage thresholds
	uint64_t total_reads_mapped_cov; // Set in 'compute_lis_alignment2'
	char filesig = '>';
	std::string suffix; // 'fasta' | 'fastq' TODO: remove?

	// Non-synchronized
	uint64_t number_total_read; // total number of reads in file. Should be known before processing and index loading. 'calculate'
	off_t    full_file_size; // the size of the full reads file (in bytes). 'calculate'
	uint64_t full_read_main; // total number of nucleotides in all reads i.e. sum of length of All read sequences 'calculate'
	std::vector<uint64_t> reads_matched_per_db; // total number of reads matched for each database
	uint64_t total_reads_denovo_clustering; // total number of reads for de novo clustering. Synchronize? - only incremented by threads.

	int otu_total; // total number of OTUs. TODO: Synchronize?, 
	std::map<string, vector<string>> otu_map; // TODO: may be find a better place for this.

	Readstats(Runopts & opts)
		:
		opts(opts),
		min_read_len(READLEN),
		max_read_len(0),
		total_reads_mapped(0),
		total_reads_mapped_cov(0),
		number_total_read(0),
		full_file_size(0),
		full_read_main(0),
		reads_matched_per_db(opts.indexfiles.size(), 0),
		total_reads_denovo_clustering(0),
		otu_total(0)
	{
		calcSuffix();
		opts.exit_early = check_file_format();
		calculate();
	}

	~Readstats() {}

	void calculate(); // calculate statistics from readsfile
	void calculate2(Runopts & opts); // old version - copy of 'compute_read_stats'
	bool check_file_format();
	void calcSuffix();
}; // ~struct Readstats
