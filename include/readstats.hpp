#pragma once
/**
 * FILE: readstats.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Collective Statistics for all Reads. Encapsulates old 'compute_read_stats' logic and results
 * Some statistics computed during Alignment, and some in Post-processing
 */

#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <atomic>

#include "common.hpp"
#include "options.hpp"

// forward
class KeyValueDatabase;

struct Readstats {
	Runopts & opts;

	// TODO: synchronize - Computed in worker thread (per read) - shared
	uint32_t min_read_len; // length of the shortest Read in the Reads file. 'parallelTraversalJob'
	uint32_t max_read_len; // length of the longest Read in the Reads file. 'parallelTraversalJob'
	std::atomic_uint64_t total_reads_mapped; // total number of reads mapped passing E-value threshold. Thread accessed: 'compute_lis_alignment'
	// thread accessed: 'compute_lis_alignment'
	std::atomic_uint64_t total_reads_mapped_cov; // total number of reads mapped passing E-value threshold and %id && %query coverage thresholds
	char filesig = FASTA_HEADER_START;
	std::string suffix; // 'fasta' | 'fastq' TODO: remove?

	// Non-synchronized
	uint64_t number_total_read; // total number of reads in file. Should be known before processing and index loading. 'calculate'
	off_t    full_file_size; // the size of the full reads file (in bytes). 'calculate'
	uint64_t full_read_main; // total number of nucleotides in all reads i.e. sum of length of All read sequences 'calculate'
	// TODO: thread accessed: 'compute_lis_alignment' - synchronize
	std::vector<uint64_t> reads_matched_per_db; // total number of reads matched for each database.
	// TODO: currently accessed in single thread ('computeStats') but potentially could be multiple threads
	uint64_t total_reads_denovo_clustering; // total number of reads for de novo clustering.

	// Clustering of reads around references by similarity i.e. 
	// {ref: [read,read,...], ref: [read,read...], ...}
	// calculated after alignment is done on all reads
	std::map<std::string, std::vector<std::string>> otu_map;

	std::mutex total_reads_mapped_cov_lock; // mutex for total_reads_mapped_cov
	std::mutex otu_map_lock; // mutex for otu_map

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
		total_reads_denovo_clustering(0)
	{
		calcSuffix();
		opts.exit_early = check_file_format();
		calculate(); // number_total_read only
	}

	~Readstats() {}

	void calculate(); // calculate statistics from readsfile
	bool check_file_format();
	void calcSuffix();
	std::string toString();
	bool restoreFromDb(KeyValueDatabase & kvdb);
	void pushOtuMap(std::string & ref_seq_str, std::string & read_seq_str);
	void printOtuMap(std::string otumapfile);
}; // ~struct Readstats
