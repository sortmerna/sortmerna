#pragma once
/**
* FILE: index.hpp
* Created: Nov 06, 2017 Mon
*/

#include <vector>

#include "options.hpp"
#include "readstats.hpp"
#include "output.hpp"
#include "indexdb.hpp"

/**
* 1. Each reference file can be indexed into multiple index parts depending on the size.
*    Each index file name follows a pattern <Name_Part> e.g. index1_0, index1_1 etc.
*/
struct Index {
	Runopts & opts;

	// Index
	uint16_t index_num = 0; // currrently loaded index number (DB file) Set in Main thread
	uint32_t part = 0; // currently loaded index part
	uint32_t number_elements = 0; /* number of positions in (L+1)-mer positions table */

	// Index stats
	long _match = 0;    /** Smith-Waterman score for a match */
	long _mismatch = 0; /** Smith-Waterman score for a mismatch */
	long _gap_open = 0; /** Smith-Waterman score for gap opening */
	long _gap_extension = 0; /** Smith-Waterman score for gap extension */

	std::vector<uint16_t> num_index_parts; /* number of parts each index file has i.e. each index can have multiple parts. See 'load_stats' */
	std::vector<vector<index_parts_stats>> index_parts_stats_vec; /** statistics for index files' parts */
	std::vector<uint64_t> full_ref;   /* corrected size of each reference index (for computing E-value) See Index::load_stats */
	std::vector<uint64_t> full_read;  /* corrected size of reads (for computing E-value) See Index::load_stats */
	std::vector<uint32_t> lnwin;      /* length of seed (sliding window L). Unique per DB. Const. Obtained in Main thread. Thread safe. See 'load_stats' */
	std::vector<uint32_t> partialwin; /* length of seed/2 */
	std::vector<uint32_t> minimal_score; /* minimal SW score in order to reach threshold E-value */
	std::vector<pair<double, double>> gumbel; // Gumbel parameters Lambda and K. Calculated in 'load_stats'
	std::vector<uint64_t> numbvs; /* number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward */
	std::vector<uint64_t> numseq;  /* total number of reference sequences in one complete reference database */

	std::vector<kmer> lookup_tbl; /**< reference to L/2-mer look up table */
	std::vector<kmer_origin> positions_tbl; /**< reference to (L+1)-mer positions table */

	Index(Runopts & opts, Readstats & readstats, Output & output)
		:
		opts(opts),
		num_index_parts(opts.indexfiles.size(), 0),
		full_ref(opts.indexfiles.size(), 0),
		full_read(opts.indexfiles.size(), readstats.full_read_main), /* ReadStats::full_read_main: total number of nucleotides in all reads <- compute_read_stats */
		lnwin(opts.indexfiles.size(), 0),
		partialwin(opts.indexfiles.size(), 0),
		minimal_score(opts.indexfiles.size(), 0),
		gumbel(opts.indexfiles.size(), std::pair<double, double>(-1.0, -1.0)),
		numbvs(opts.indexfiles.size(), 0),
		numseq(opts.indexfiles.size(), 0)
	{
		load_stats(readstats, output);
	}
	~Index() {}

	// args: index number and number of the part of the index
	void load(uint32_t idx_num, uint32_t idx_part);
	void load_stats(Readstats & readstats, Output & output); // TODO: make private?
}; // ~struct Index
