/*
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * Copyright (C) 2014 Bonsai Bioinformatics Research Group
 *
 * This file is part of SortMeRNA.
 *
 * SortMeRNA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SortMeRNA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * file: paralleltraversal.hpp
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */

#ifndef PARALLELTRAVERSAL_H
#define PARALLELTRAVERSAL_H

#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>
#include <queue>
#include "common.hpp"
#include "bitvector.hpp"
#include "indexdb.hpp"
#include <errno.h>



using namespace std;


extern timeval t;

/* for each 18-mer hit on the read, we store the 
   key to find the positions and the window number
   on the read at which the 18-mer occurs   */
struct id_win
{
	/// key value to find index positions
	uint32_t id;
	/// the associated window number on the read 
	uint32_t win;
};


typedef pair<uint32_t,uint32_t> mypair;



	
/*
 *
 * FUNCTION	: paralleltraversal 
 *
 * PARAMETERS	: char* inputreads			- reads input file
 *		  		  	char* ptr_filetype_ar	- accepted reads output file
 *		  		  	char* ptr_filetype_or	- rejected reads output file
 *		  		  	long int match				- reward for a match (Smith-Waterman) (must be positive)
 *				  		long int mismatch			- penalty for a mismatch (Smith-Waterman) (must be negative)
 *				  		long int gap_open			- penalty for a gap (Smith-Waterman) (must be positive, will be used as negative)
 *							long int gap_extenstion	- penalty for extending a gap (Smith-Waterman) (must be positive, will be used as negative)
 *							string &root					- reference database filename without extension (prefix used in indexdb)
 *							string &path					- path to sortmedna index directory
 *							int argc							- number of arguments in command used to launch sortmedna (for SAM file)
 *							char** argv						- command used to launch sortmedna (for SAM file)
 *
 * OUTPUT	 : none
 * 
 **************************************************************************************************************/
void
paralleltraversal ( char* inputreads,
                   char* ptr_filetype_ar,
                   char* ptr_filetype_or,
                   int32_t match,
                   int32_t mismatch,
                   int32_t gap_open,
                   int32_t gap_extension,
                   int32_t score_N,
                   vector< vector<uint32_t> >& skiplengths,
                   int argc,
                   char **argv,
                   bool yes_SQ,
                   vector< pair<string,string> >& myfiles,
                   bool exit_early);


/*
 *
 * FUNCTION	: traversetrie_align()
 *		  exact matching of [p_1] in [w_1] is completed fully in the trie nodes, continue parallel 
 *		  traversal of BTRIE beginning at [w_2]
 *                  
 *		  	w = |------ [w_1] ------|------ [w_2] ------|
 *		  	p = |------ [p_1] ------|------ [p_2] --....--|	 
 *		  	    |------ trie ------x|----- tail ----....--| 
 *
 * PARAMETERS   : NodeElement* curr		- the last trie node 'x' of exact matching 
 *		  unsigned short int root_l	- initial levenshtein state, root_l = 0
 *		  unsigned char depth			- the depth 'x' in the BTRIE
 *
 * OUTPUT	: none
 *
 * NOTES	: final accepting states for k = 1 are [8,13]
 *					     k = 2 are [50,89]
 *					     k = 3 are [322,601]
 */

inline void
traversetrie_align ( NodeElement *trie_t,
                    uint32_t lev_t,
                    unsigned char depth,
                    MYBITSET *win_k1_ptr,
                    MYBITSET *win_k1_full,
                    bool &accept_zero_kmer,
                    vector< id_win > &id_hits,
                    uint32_t readn,
                    uint32_t win_num,
                    uint32_t partialwin
                    );


void 
traverse_btrie ( NodeElement* trie_node, unsigned int &hash_id, unsigned char depth );


void preprocess_data(vector< pair<string,string> >& myfiles,
                     char** argv,
                     int argc,
                     bool yes_SQ,
                     char* acceptedstrings_sam,
                     int32_t _match,
                     int32_t _mismatch,
                     int32_t _gap_open,
                     int32_t _gap_extension,
                     vector<vector<uint32_t> >& skiplengths,
                     vector<uint16_t>& num_index_parts,
                     vector<vector<index_parts_stats> >& index_parts_stats_vec,
                     vector<uint64_t>& full_ref,
                     vector<uint64_t>& full_read,
                     vector<uint32_t>& lnwin,
                     vector<uint32_t>& partialwin,
                     vector<uint32_t>& minimal_score,
                     uint32_t number_total_read,
                     vector<pair<double, double> >& gumbel,
                     vector<uint32_t>& numbvs,
                     vector<uint32_t>& numseq);

void
load_index( char* ptr_dbindex,
           string part_str,
           kmer*& lookup_tbl,
           kmer_origin*& positions_tbl,
           uint32_t& number_elements,
           uint32_t lnwin);

void
find_lis( deque<mypair> &a, vector<int> &b );


#endif //~parallel_traversal_h
