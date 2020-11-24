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
 * 2. 'total_mapped_sw_id_cov'
 *        Calculated during alignment and stored to KVDB
 *        Thread accessed - Synchronize
 * 3. 'reads_matched_per_db' - Synchronize.
 *			Calculated during alignment. Thread accessed.
 * 4. 'total_reads_denovo_clustering'
 *			TODO: currently accessed in single thread ('computeStats') but potentially could be multiple threads
 */
struct Readstats 
{
	std::string dbkey; // Hashed concatenation of underscore separated basenames of the read files. Used as the key into the Key-value DB. 
	std::string suffix; // 'fasta' | 'fastq' TODO: remove?

	std::atomic<uint64_t> total_aligned; // reads passing E-value threshold.
	std::atomic<uint64_t> total_aligned_id_cov; // [2] reads passing E-value, %ID, %COV (query coverage) thresholds.
	std::atomic<uint64_t> total_denovo; // [4] 'de novo' reads (SW + !ID + !COV).
	std::atomic<uint64_t> short_reads_num; // reads shorter than a threshold of N nucleotides. Reset for each index.

	uint64_t all_reads_count; // [1] total number of reads in file.
	uint64_t all_reads_len; // total number of nucleotides in all reads i.e. sum of length of All read sequences
	uint64_t total_otu;
    uint32_t min_read_len; // shortest Read length. (read only)
    uint32_t max_read_len; // longest Read length. (read only)

	std::vector<uint64_t> reads_matched_per_db; // [3] reads matched per database.
    //              |_TODO: should be atomic std::atomic<uint64_t> 20201019

	bool is_stats_calc; // flags 'computeStats' was called.
	bool is_set_aligned_id_cov; // flag 'total_aligned_id_cov' was calculated (so no need to calculate no more)

	Readstats(uint64_t all_reads_count, uint64_t all_reads_len, uint32_t min_read_len, uint32_t max_read_len, KeyValueDatabase& kvdb, Runopts& opts);

	void calcSuffix(Runopts& opts);
	std::string toBstring();
	std::string toString();
	bool restoreFromDb(KeyValueDatabase & kvdb);
	void store_to_db(KeyValueDatabase & kvdb);
	void set_is_set_aligned_id_cov();
}; // ~struct Readstats
