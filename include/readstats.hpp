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

/*
 * 1. 'all_reads_count' - Should be known before processing and index loading. 
 * 2. 'total_reads_mapped_cov' Synchronize. Thread accessed from 'compute_lis_alignment'
 * 3. 'reads_matched_per_db' - Synchronize.
 *			Calculated in 'compute_lis_alignment' during alignment. Thread accessed.
 * 4. 'total_reads_denovo_clustering'
 *			TODO: currently accessed in single thread ('computeStats') but potentially could be multiple threads
 *          Setter: 'computeStats'.
 *			User: 'writeLog'
 * 5. 'otu_map' - Clustering of reads around references by similarity i.e. {ref: [read,read,...], ref: [read,read...], ...}
 *			calculated after alignment is done on all reads
 *			Setter: 'computeStats' post-processing callback
 *			User: 'printOtuMap'
 *			TODO: Store in DB? Can be very big.
 */
struct Readstats 
{
	std::string dbkey; // key to store/retrieve reads stats in Key-value DB. Set once at construct time.
	char filesig = FASTA_HEADER_START;
	std::string suffix; // 'fasta' | 'fastq' TODO: remove?

	std::atomic<uint32_t> min_read_len; // length of the shortest Read in the Reads file. 'parallelTraversalJob'
	std::atomic<uint32_t> max_read_len; // length of the longest Read in the Reads file. 'parallelTraversalJob'
	std::atomic<uint64_t> total_reads_mapped; // total number of reads mapped passing E-value threshold. Set in 'compute_lis_alignment'
	std::atomic<uint64_t> total_reads_mapped_cov; // [2] total number of reads mapped passing E-value, %id, %query coverage thresholds

	uint64_t all_reads_count; // [1] total number of reads in file. Non-sync. 'Readstats::calculate'
	uint64_t all_reads_len; // total number of nucleotides in all reads i.e. sum of length of All read sequences 'Readstats::calculate'
	uint64_t total_reads_denovo_clustering; // [4] total number of reads for de novo clustering. 'computeStats' post-processing callback

	std::vector<uint64_t> reads_matched_per_db; // [3] total number of reads matched for each database. `compute_lis_alignment`.
	std::map<std::string, std::vector<std::string>> otu_map; // [5] Populated in 'computeStats' post-processor callback

	bool stats_calc_done; // flags 'computeStats' was called. Set in 'postProcess'

	Readstats(Runopts & opts, KeyValueDatabase &kvdb);
	~Readstats() {}

	void calculate(Runopts &opts); // calculate statistics from readsfile
	void calcSuffix(Runopts &opts);
	std::string toBstring();
	std::string toString();
	bool restoreFromDb(KeyValueDatabase & kvdb);
	void store_to_db(KeyValueDatabase & kvdb);
	void pushOtuMap(std::string & ref_seq_str, std::string & read_seq_str);
	void printOtuMap(std::string otumapfile);
}; // ~struct Readstats
