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
 * along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 * @endparblock
 *
 * @authors jenya.kopylov@gmail.com
 *          laurent.noe@lifl.fr
 *          helene.touzet@lifl.fr
 *          robknight@ucsd.edu
 */

 /** @file paralleltraversal.hpp */

#ifndef PARALLELTRAVERSAL_H
#define PARALLELTRAVERSAL_H

#include "outputformats.hpp"
#include "load_index.hpp"
#include "traverse_bursttrie.hpp"

#include <iomanip>
#include <map>
#include <algorithm>
#include <queue>
 

using namespace std;

/*! @brief Number of slots by which to dynamically
           increment the array storing all alignments
           per read */
#define BEST_HITS_INCREMENT 100

/* holds the index of the minimum and maximum scoring
   alignments in an array of alignments pointed to by
   s_align* ptr */
struct alignment_struct
{
  uint32_t max_size;
  uint32_t size;
  uint32_t min_index;
  uint32_t max_index;
  s_align* ptr;
  alignment_struct(uint32_t max_size, uint32_t size, uint32_t min, uint32_t max, s_align* p) : max_size(max_size), size(size), min_index(min), max_index(max), ptr(p) {}
};


/*! @brief Type mypair

    A data structure holding two variables
    of type uint32_t.
*/
typedef pair<uint32_t,uint32_t> mypair;

/*! @fn paralleltraversal()
    @brief Traverse the query input and indexed database and output
           alignments passing the E-value threshold
    @detail The following methods will be executed:
    <ol> 
      <li> divide large read files into mmap'd regions,
           taking into account the read (and its pair) which may
           be split between two file sections </li>
      <li> compute the gumbel parameters (lamda and K) using ALP,
           load the index fully or in parts (depending on how
           it was built) </li>
      <li> using 3 intervals, scan over the read and collect all
           L-mers on the read which match to the reference
           index with at most 1 error. This is done using
           parallel traversal between the index and the
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

    @param char* inputreads, pointer to query reads file
    @param *ptr_filetype_ar, pointer to string for aligned seqeunces filepath
    @param *ptr_filetype_or, pointer to string for rejected sequences filepath
    @param int32_t match, SW match reward score (positive)
    @param int32_t mismatch, SW mismatch penalty score (negative)
    @param int32_t gap_open, SW gap open penalty score (positive)
    @param int32_t gap_extension, SW gap extend penalty score (positive)
    @param int32_t score_N, SW penalty for ambiguous nucleotide (negative)
    @param vector< vector<uint32_t> >& skiplengths, 
    @param int argc, number of arguments passed to sortmerna
    @param char **argv, argument string passed to sortmerna
    @param bool yes_SQ, boolean to include @SQ tags in SAM output
    @param vector< pair<string,string> >& myfiles, 
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
                    bool exit_early );

void
find_lis( deque<mypair> &a, vector<int> &b );


#endif //~parallel_traversal_h
