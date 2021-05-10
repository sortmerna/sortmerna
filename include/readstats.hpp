/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

/*
 * file: readstats.hpp
 * created: Nov 06, 2017 Mon
 *
 * Collective Statistics for all Reads. Encapsulates old 'compute_read_stats' logic and results
 * Some statistics computed during Alignment, and some in Post-processing
 */

#pragma once

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

	uint64_t all_reads_count; // [1] total number of reads in file.
	uint64_t all_reads_len; // total number of nucleotides in all reads i.e. sum of length of All read sequences
	uint32_t min_read_len; // shortest Read length. (read only)
	uint32_t max_read_len; // longest Read length. (read only)
	uint64_t total_otu; // not to store in DB

	std::atomic<uint64_t> num_aligned; // reads passing E-value threshold.
	std::atomic<uint64_t> n_yid_ncov; // SW + ID - COV i.e. aligned passing ID, failing COV
	std::atomic<uint64_t> n_nid_ycov; // SW - ID + COV i.e. aligned failing ID, passing COV
	std::atomic<uint64_t> n_yid_ycov; // [2] SW + ID + COV i.e. aligned passing ID, passing COV
	std::atomic<uint64_t> num_denovo; // [4] SW - ID - COV i.e. 'de novo' reads, aligned failing ID, failing COV
	std::atomic<uint64_t> num_short; // count of reads shorter than a threshold of N nucleotides. Reset for each index.

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
