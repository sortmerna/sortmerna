#pragma once
/**
* FILE: index.hpp
* Created: Nov 06, 2017 Mon
*/

#include <vector>
#include <cstdint>

// forward
struct Runopts;
struct kmer;
struct kmer_origin;
class Refstats;

/**
 * 1. Each reference file can be indexed into multiple index parts depending on the file size.
 *    Each index file name follows a pattern <Name_Part> e.g. index1_0, index1_1 etc.
 */
struct Index {
	uint16_t index_num = 0; // currrently loaded index number (DB file) Set in Main thread
	uint32_t part = 0; // currently loaded index part
	uint32_t number_elements = 0; /* number of positions in (L+1)-mer positions table */

	std::vector<kmer> lookup_tbl; /**< reference to L/2-mer look up table */
	std::vector<kmer_origin> positions_tbl; /**< reference to (L+1)-mer positions table */

	// Index stats
	//long _match = 0;    /* Smith-Waterman score for a match */
	//long _mismatch = 0; /* Smith-Waterman score for a mismatch */
	//long _gap_open = 0; /* Smith-Waterman score for gap opening */
	//long _gap_extension = 0; /* Smith-Waterman score for gap extension */

	Index(Runopts & opts);
	~Index() {}

	void load(uint32_t idx_num, uint32_t idx_part, Runopts & opts, Refstats & refstats);
	void clear();
}; // ~struct Index
