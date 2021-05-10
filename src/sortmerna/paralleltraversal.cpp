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

/**
 * @file paralleltraversal.cpp
 * @brief functions for index traversal.
 */

#include <algorithm>
#include <locale>
#include <iomanip> // output formatting

#include "paralleltraversal.hpp"
#include "kseq.h"
#include "kseq_load.hpp"
#include "traverse_bursttrie.hpp"
#include "alignment.hpp"

#include "options.hpp"
#include "ThreadPool.hpp"
#include "read.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "index.hpp"
#include "references.hpp"
#include "kvdb.hpp"
#include "processor.hpp"
#include "readfeed.hpp"
#include "output.hpp"
#include "readsqueue.hpp"


#if defined(_WIN32)
#define O_SMR_READ_BIN O_RDONLY | O_BINARY
#else
#define O_SMR_READ_BIN O_RDONLY
#endif

// forward
int clear_dir(std::string dpath);

 // see "heuristic 1" below
 //#define HEURISTIC1_OFF


/* 
 * Callback run in a Processor thread
 * Called on each index * index_part * read.num_strands
 *
 * @param isLastStrand flags when the last strand (out of max 2 strands) is passed for matching
 */
void traverse
	(
		Runopts& opts, 
		Index& index, 
		References& refs, 
		Readstats& readstats, 
		Refstats& refstats, 
		Read& read,
		bool isLastStrand
	)
{
	read.lastIndex = index.index_num;
	read.lastPart = index.part;

	bool read_to_count = true; // passed directly to compute_lis_alignment. TODO: What's the point?

	uint32_t win_shift = opts.skiplengths[index.index_num][0];
	// keep track of windows (read positions) which have already been traversed
	// in the burst trie using different shifts. Initially all False
	vector<bool> read_pos_searched(read.sequence.size());

	size_t pass_n = 0; // Pass number (possible value 0,1,2)
	uint32_t max_SW_score = read.sequence.size() * opts.match; // the maximum SW score attainable for this read

	std::vector<UCHAR> bitvec; // window (prefix/suffix) bitvector

	// TODO: below 2 values are unique per index part. Move to index?
	uint32_t bitvec_size = (refstats.partialwin[index.index_num] - 2) << 2; // e.g. 9 - 2 = 0000 0111 << 2 = 0001 1100 = 28
	// Does this mark where in 32-bit the bitvector starts?
	uint32_t offset = (refstats.partialwin[index.index_num] - 3) << 2; // e.g. 9 - 3 = 0000 0110 << 2 = 0001 1000 = 24

	// loop search positions on the read in multiple passes
	// changing the step (skip length/windowshift) when necessary
	for (bool search = true; search; )
	{
		// number of k-mer windows fit along the read given 
		// the window size and the search step (windowshift)
		uint32_t numwin = ( 
				read.sequence.size() - refstats.lnwin[index.index_num] + win_shift
			) / win_shift;

		uint32_t win_pos = 0; // position (index) of the window's first char on the read i.e. [0...read.sequence.length-1]
		// iterate the windows
		for (uint32_t win_num = 0; win_num < numwin; ++win_num)
		{
			if (read.is04) read.flip34(); // Make sure the read is in 03 encoding for index search

			// skip position when the seed at this position has already been searched for in a previous Passes
			if (!read_pos_searched[win_pos])
			{
				read_pos_searched[win_pos] = true; // mark position as searched
				// this flag it set to true if a match is found during
				// subsearch 1(a), to skip subsearch 1(b)
				bool accept_zero_kmer = false;
				// ids for k-mers hits on the reference database
				vector<id_win> id_hits; // TODO: add directly to 'id_win_hits'? - No, id_win_hits may contain hits from different index parts.

				bitvec.resize(bitvec_size);
				std::fill(bitvec.begin(), bitvec.end(), 0);

				auto ii = win_pos + refstats.partialwin[index.index_num];
				init_win_f(&read.isequence[ii],	&bitvec[0],	&bitvec[4],	refstats.numbvs[index.index_num]);

				// the hash of the 'first half' of the kmer window
				uint32_t keyf = read.hashKmer(win_pos, refstats.partialwin[index.index_num]);

				// TODO: remove in production
				if (index.lookup_tbl.size() <= keyf) {
					size_t vsize = index.lookup_tbl.size();
					uint16_t idxn = index.index_num;
					uint16_t idxp = index.part;
					std::string id = read.id;
					bool is03 = read.is03;
					bool is04 = read.is04;
					ERR("lookup index: ", keyf, " is larger than lookup_tbl.size: ", vsize, 
						" Index: ", idxn, " Part: ", idxp, " Read.id: ", id, " Read.is03: ", is03, " Read.is04: ", is04, " Aborting..");
					exit(EXIT_FAILURE);
				}

				// do traversal if the exact half window exists in the burst trie
				if ( index.lookup_tbl[keyf].count > opts.minoccur && index.lookup_tbl[keyf].trie_F != NULL )
				{
					/* subsearch (1)(a) d([p_1],[w_1]) = 0 and d([p_2],[w_2]) <= 1;
					*
					*  w = |------ [w_1] ------|------ [w_2] ------|
					*  p = |------ [p_1] ------|------ [p_2] ----| (0/1 deletion in [p_2])
					*              or
					*    = |------ [p_1] ------|------ [p_2] ------| (0/1 match/substitution in [p_2])
					*        or
					*    = |------ [p_1] ------|------ [p_2] --------| (0/1 insertion in [p_2])
					*
					*/
					traversetrie_align(
						index.lookup_tbl[keyf].trie_F,
						0,
						0,
						&bitvec[0],
						&bitvec[offset],
						accept_zero_kmer,
						id_hits,
						win_pos,
						refstats.partialwin[index.index_num],
						opts
					);
				} //~if exact half window exists in the burst trie

				// only search rear kmer if an exact match has not been found for the forward
				if (!accept_zero_kmer)
				{
					//bitvec.resize(bitvec_size);
					std::fill(bitvec.begin(), bitvec.end(), 0);

					// init the first bitvector window
					auto ii = win_pos + refstats.partialwin[index.index_num] - 1;
					init_win_r(&read.isequence[ii],	&bitvec[0],	&bitvec[4],	refstats.numbvs[index.index_num]);

					// the hash of the second (rear) half of the kmer window
					uint32_t keyr = read.hashKmer(win_pos + refstats.partialwin[index.index_num], refstats.partialwin[index.index_num]);

					// TODO: remove in production
					if (index.lookup_tbl.size() <= keyr) {
						size_t vsize = index.lookup_tbl.size();
						uint16_t idxn = index.index_num;
						uint16_t idxp = index.part;
						std::string id = read.id;
						bool is03 = read.is03;
						bool is04 = read.is04;
						ERR("Thread: ", std::this_thread::get_id(), " lookup index: ", keyr, 
							" is larger than lookup_tbl.size: ", vsize, " Index: ", idxn, " Part: ", idxp, 
							" Read.id: ", id, " Read.is03: ", is03, " Read.is04: ", is04, " Aborting...");
						exit(EXIT_FAILURE);
					}

					// continue subsearch (1)(b)
					if ( index.lookup_tbl[keyr].count > opts.minoccur && index.lookup_tbl[keyr].trie_R != NULL )
					{
						/* subsearch (1)(b) d([p_1],[w_1]) = 1 and d([p_2],[w_2]) = 0;
						*
						*  w =    |------ [w_1] ------|------ [w_2] -------|
						*  p =      |------- [p_1] ---|--------- [p_2] ----| (1 deletion in [p_1])
						*              or
						*    =    |------ [p_1] ------|------ [p_2] -------| (1 match/substitution in [p_1])
						*        or
						*    = |------- [p_1] --------|---- [p_2] ---------| (1 insertion in [p_1])
						*
						*/
						traversetrie_align(
							index.lookup_tbl[keyr].trie_R,
							0,
							0,
							&bitvec[0],
							&bitvec[offset],
							accept_zero_kmer,
							id_hits,
							win_pos,
							refstats.partialwin[index.index_num], 
							opts);
					}//~if exact half window exists in the reverse burst trie                    
				}//~if (!accept_zero_kmer)

				// store found seed hits in the read
				if (!id_hits.empty())
				{
					for (auto const& hit: id_hits)
					{
						read.id_win_hits.push_back(hit);
					}
					++read.hit_seeds;
				}
			} // ~if not read_pos_searched[win_pos]

			// all k-mers for a given shift-size are to be looked up prior proceeding to the LIS/SW calculation
			if (win_num == numwin - 1)
			{
				// calculate LIS if the number of matching seeds on the read meets the threshold (default 2)
				if (read.hit_seeds >= (uint32_t)opts.num_seeds) {
					compute_lis_alignment(read, opts, index, refs, readstats, refstats,	search,	max_SW_score);
				}

				// if the read was not accepted at the current shift,
				// use the next (smaller) window shift
				if (search)
				{
					if (pass_n == 2) 
						search = false; // the last (3rd) Pass has been made
					else
					{
						// the next interval size equals to the current one, skip it
						while (pass_n < 3
							&& opts.skiplengths[index.index_num][pass_n] == opts.skiplengths[index.index_num][pass_n + 1])
							++pass_n;
						if (++pass_n > 2) search = false;
						// set interval skip length for next Pass
						else win_shift = opts.skiplengths[index.index_num][pass_n];
					}
				}
				break; // go to the next shift size
			}//~( win_num == NUMWIN-1 )
			win_pos += win_shift;
		}//~for (each window)                
			//~while all skip/shift lengths have not been tested, or a match has not been found
	}// ~while (search);

	// all_N_best_max_SW Or all_N hits found - stop further processing of this read
	if (opts.num_alignments > 0) {
		if ((opts.is_best && opts.num_alignments == read.max_SW_count) ||
			(!opts.is_best && read.alignment.alignv.size() == opts.num_alignments)) {
			read.is_done = true;
		}
	}
	// end of processing and read.alignments > 0
	else {
		bool is_last_idx = (index.index_num == opts.indexfiles.size() - 1) && (index.part == refstats.num_index_parts[index.index_num] - 1);
		if (is_last_idx && isLastStrand && read.alignment.alignv.size() > 0)
			read.is_done = true;
	}
} // ~traverse

/**
 * verify the alignment was already performed by querying the KVDB
 * Alignment descriptor:
 *   List all index files: hash, size
 *   List read files: hash, size
 *   List reference files: hash, size
 *   List number of aligned reads
 *   Store options and compare to the current. Add '==' operator.
 *   Store the list of DBKeys of all aligned reads (?)
 *
 * Alignment IS Done IF
 *  - reads files are the same
 *  - references are the same
 *  - index is present and the names/hashes are the same as stored in DB
 *  - alignment results are stored
 *  - read statistics are stored and is_done = True
 */
bool is_aligned(Runopts& opts, Readstats& readstats, Output& output, Index& index, KeyValueDatabase& kvdb)
{
	INFO("TODO");
	return false;
}