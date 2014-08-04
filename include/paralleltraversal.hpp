/**
 * @file paralleltraversal.hpp
 * @brief Function and variable definitions for paralleltraversal.cpp
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright Copyright (C) 2012-2014 Bonsai Bioinformatics Research Group, LIFL and 
 * INRIA Nord-Europe, France
 * OTU-picking extensions developed in the Knight Lab, BioFrontiers Institute,
 * University of Colorado at Boulder, Boulder, CO
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
 * @endparblock
 *
 * @authors jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 */

 /** @file paralleltraversal.hpp */

#ifndef PARALLELTRAVERSAL_H
#define PARALLELTRAVERSAL_H

#include "bitvector.hpp"
#include "outputformats.hpp"
//! ALP program for computing the Gumbel parameters
#include "../alp/sls_alp_data.hpp"
#include "../alp/sls_alp_sim.hpp"
#include "../alp/gumbel_params.hpp"

#include <iomanip>
#include <map>
#include <algorithm>
#include <queue>
#include <deque> 


using namespace std;


//extern timeval t;

/* for each 18-mer hit on the read, we store the 
   key to find the positions and the window number
   on the read at which the 18-mer occurs   */
struct id_win
{
	// key value to find index positions
	uint32_t id;
	// the associated window number on the read 
	uint32_t win;
};


/*! @brief Type mypair

    A data structure holding two variables
    of type uint32_t.
*/
typedef pair<uint32_t,uint32_t> mypair;



	
/*! @fn paralleltraversal()
    @brief traverse the query input and indexed database and output
           alignments passing the E-value threshold
    @detail The main function of SortMeRNA with the following methods:
    <ol> 
      <li> divide large read files into mmap'd regions,
           taking into account the read (and its pair) which may
           be split between two file sections </li>
      <li> load the index, compute the gumbel parameters
           (lamda and K) using ALP </li>
      <li> using 3 intervals, scan over the read and collect all
           L-mers on the read which match to the
           reference index with at most 1 error. This is done 
           using parallel traversal between the index and the
           Levenshtein automaton </li>
      <li> if enough L-mers were collected, extend them into
           longer matches using the Longest Increasing
           subsequence (LIS) of positions where the L-mers
           matched on the read and the associated reference
           sequences </li>
      <li> if the LIS is long enough, use the starting positions
           of the LIS to estimate the starting position
           of an alignment and pass this reference segment and
           read to SSW </li>
      <li> if the alignment score is at least the minimum score
           corresponding to the E-value threshold, keep the read,
           otherwise continue searching for other LIS or more
           L-mers using smaller intervals </li>
    </ol>

    @param char* inputreads
    @param *ptr_filetype_ar
    @param *ptr_filetype_or
    @param int32_t match
    @param int32_t mismatch
    @param int32_t gap_open
    @param int32_t gap_extension
    @param int32_t score_N
    @param vector< vector<uint32_t> >& skiplengths
    @param int argc
    @param char **argv
    @param bool yes_SQ
    @param vector< pair<string,string> >& myfiles
    @return void
    @version 1.0 Jan 14, 2013 
*/
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
