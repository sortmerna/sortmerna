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
 * file: index.hpp
 * created: Nov 06, 2017 Mon
 */

#pragma once

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
	uint16_t index_num; // currrently loaded index number (DB file) Set in Main thread
	uint32_t part; // currently loaded index part
	uint32_t number_elements; /* number of positions in (L+1)-mer positions table */
	bool is_ready; // flags the index is built and ready

	// Index stats
	//long _match = 0;    /* Smith-Waterman score for a match */
	//long _mismatch = 0; /* Smith-Waterman score for a mismatch */
	//long _gap_open = 0; /* Smith-Waterman score for gap opening */
	//long _gap_extension = 0; /* Smith-Waterman score for gap extension */

	std::vector<kmer> lookup_tbl; /**< reference to L/2-mer look up table */
	std::vector<kmer_origin> positions_tbl; /**< reference to (L+1)-mer positions table */

	/*
	 * Initilize the index.
	 * If index files do not exist or are empty - build the index.
	 */
	Index(Runopts & opts);
	//~Index() {}
	void load(uint32_t idx_num, uint32_t idx_part, std::vector<std::pair<std::string, std::string>>& indexfiles, Refstats & refstats);
	void unload();
}; // ~struct Index
