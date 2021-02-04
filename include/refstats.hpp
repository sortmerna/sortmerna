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
 * file: refstats.hpp
 * created: Dec 23, 2017 Sat
 */

#pragma once

#include <cstdint>
#include <vector>
#include <utility> // std::pair

#include "indexdb.hpp" // index_parts_stats;

// forward
struct Readstats;
struct Runopts;


class Refstats {
public:
	std::vector<uint16_t> num_index_parts; /* number of parts in each index file (index can have multiple parts). see Refstats::load */
	std::vector<std::vector<index_parts_stats>> index_parts_stats_vec; /* index parts statistics */
	std::vector<uint64_t> full_ref;   /* corrected size of each reference index (for computing E-value) see Refstats::load */
	std::vector<uint64_t> full_read;  /* corrected size of reads (for computing E-value) see Refstats::load */

	/* 
	 * length of seed (sliding window L). Unique per DB. Const. 
	 * Set by 'opts.seed_win_len', used in indexing, stored in '.stats' table, and loaded here by 'load'
	 */
	std::vector<uint32_t> lnwin;

	/* length of seed/2 */
	std::vector<uint32_t> partialwin;

	/* minimal SW score corresponing to the threshold E-value */
	std::vector<uint32_t> minimal_score;
	std::vector<std::pair<double, double>> gumbel; // Gumbel parameters Lambda and K. see Refstats::load
	std::vector<uint64_t> numbvs; /* number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward */
	std::vector<uint64_t> numseq;  /* total number of reference sequences in one complete reference database */

public:
	Refstats(Runopts& opts, Readstats& readstats);
	//~Refstats() {}

private:
	void load(Runopts& opts, Readstats& readstats); // called at construction
};
