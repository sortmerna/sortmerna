#pragma once
/*
 * FILE: refstats.hpp
 * Created: Dec 23, 2017 Sat
 */

#include <cstdint>
#include <vector>
#include <utility> // std::pair

#include "indexdb.hpp" // index_parts_stats;

// forward
struct Readstats;
struct Runopts;


class Refstats {
public:
	std::vector<uint16_t> num_index_parts; /* number of parts in each index file i.e. each index can have multiple parts. <--load */
	std::vector<std::vector<index_parts_stats>> index_parts_stats_vec; /* index parts statistics */
	std::vector<uint64_t> full_ref;   /* corrected size of each reference index (for computing E-value) <--load */
	std::vector<uint64_t> full_read;  /* corrected size of reads (for computing E-value) <--load */
	std::vector<uint32_t> lnwin;      /* length of seed (sliding window L). Unique per DB. Const. Obtained in Main thread. Thread safe. See 'load_stats' */
	std::vector<uint32_t> partialwin; /* length of seed/2 */
	std::vector<uint32_t> minimal_score; /* minimal SW score in order to reach threshold E-value */
	std::vector<std::pair<double, double>> gumbel; // Gumbel parameters Lambda and K. <--'load_stats'
	std::vector<uint64_t> numbvs; /* number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward */
	std::vector<uint64_t> numseq;  /* total number of reference sequences in one complete reference database */

public:
	Refstats(Runopts & opts, Readstats & readstats);
	~Refstats() {}

private:
	void load(Runopts & opts, Readstats & readstats); // called at constructions
};
