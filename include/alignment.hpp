#pragma once
/**
 * @FILE: alignment.hpp
 * @brief Function and variable definitions for alignment.cpp
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2012-16 Bonsai Bioinformatics Research Group
 * @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
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
 * @contributors Jenya Kopylova   jenya.kopylov@gmail.com
 *               Laurent Noé      laurent.noe@lifl.fr
 *               Pierre Pericard  pierre.pericard@lifl.fr
 *               Daniel McDonald  wasade@gmail.com
 *               Mikaël Salson    mikael.salson@lifl.fr
 *               Hélène Touzet    helene.touzet@lifl.fr
 *               Rob Knight       robknight@ucsd.edu
 */

#include <map>
#include <queue>
#include <algorithm>

#include "traverse_bursttrie.hpp"
#include "ssw.hpp"

 // forward
class Read;
struct Runopts;
struct Index;
class References;
class Output;
struct Readstats;
class Refstats;

using namespace std;

/*! @brief Number of slots by which to dynamically
           increment the array storing all alignments
           per read */
#define BEST_HITS_INCREMENT 100

/*! @brief Euler's constant */
#define EXP 2.71828182845904523536

/*! @brief Type mypair

    A data structure holding two variables
    of type uint32_t.
*/
typedef pair<uint32_t,uint32_t> uint32pair;


/*! @fn smallest()
    @brief Determine the smallest integer.
    @details The mypair data structure holds two integers,
             the first being the position a k-mer occurs 
             on the reference sequence and the second
             being the position a k-mer occurs on the
             query sequence. This function takes two
             mypair data structures and returns the 

    @param const pair<uint32_t,uint32_t> &a
    @param const pair<uint32_t,uint32_t> &b
    @return smallest integer of a and b, or a if a == b
*/
//inline bool smallest ( const uint32pair &a, const uint32pair &b );

/*! @fn largest()
    @brief Return the largest integer of two input integers
    @param const mypair &a
    @param const mypair &b
    @return 'a' goes before 'b' if a.first > b.first, otherwise
            'a'
*/
//bool largest ( const uint32pair &a, const uint32pair &b );

/*! @fn find_lis()
 *  @brief Given a list of matching positions on the read, find the longest
           strictly increasing subsequence, O(n log k)
    @param deque<pair<uint32_t, uint32_t> > &a  list of matching positions on the read which fall within a range of the read's length on the genome
    @param vector<uint32_t> &b  array of starting positions of each longest subsequence
*/
void find_lis(deque<pair<uint32_t, uint32_t> > &a, vector<uint32_t> &b);

/*! @brief struct alignment_struct
   holds the index of the minimum and maximum scoring
   alignments in an array of alignments pointed to by
   s_align* ptr */
struct alignment_struct
{
  uint32_t max_size; // max size of s_align array
  uint32_t size; // actual size of s_align array
  uint32_t min_index;
  uint32_t max_index;
  s_align* ptr;
  alignment_struct(): max_size(0), size(0), min_index(0), max_index(0), ptr(0) {}
  alignment_struct(uint32_t max_size,
                   uint32_t size,
                   uint32_t min,
                   uint32_t max,
                   s_align* p) : max_size(max_size), size(size), min_index(min), max_index(max), ptr(p) {}
};
/*
 * called on each idx * part * read * strand * [1..max opts.skiplengths[index_num].size (3 by default)]
 *
 * @param search OUT boolean
 *        return 'True' to indicate keep searching for more seed matches and better alignment.
 *		  return 'False' - stop search, the alignment is found
 * @param max_SW_score  the maximum SW score attainable for this read i.e. perfect match
 * @param read_to_count
 */
void compute_lis_alignment(
	Read & read, Runopts & opts, Index & index, References & refs, Readstats & readstats, Refstats & refstats,
	bool & search,
	uint32_t max_SW_score,
	bool& read_to_count
);