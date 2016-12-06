/**
 * @file alignment.hpp
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
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 */

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <map>
#include <queue>
#include <algorithm>
#include "traverse_bursttrie.hpp"
#include "outputformats.hpp"

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
typedef pair<uint32_t,uint32_t> mypair;


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
inline bool
smallest ( const mypair &a, const mypair &b );

/*! @fn largest()
    @brief Return the largest integer of two input integers
    @param const mypair &a
    @param const mypair &b
    @return 'a' goes before 'b' if a.first > b.first, otherwise
            'a'
*/
bool
largest ( const mypair &a, const mypair &b );

/*! @fn find_lis()
 *  @brief Given a list of matching positions on the read, find the longest
           strictly increasing subsequence, O(n log k)
    @param deque<pair<uint32_t, uint32_t> > &a
    @param vector<uint32_t> &b
    @param int64_t readn
*/
void
find_lis(deque<pair<uint32_t, uint32_t> > &a /**< list of matching positions on the read which fall within a range of the read's length on the genome */,
         vector<uint32_t> &b /**< array of starting positions of each longest subsequence */,
         uint64_t readn /**< read number */);

/*! @brief struct alignment_struct
   holds the index of the minimum and maximum scoring
   alignments in an array of alignments pointed to by
   s_align* ptr */
struct alignment_struct
{
  uint32_t max_size;
  uint32_t size;
  uint32_t min_index;
  uint32_t max_index;
  s_align* ptr;
  alignment_struct(uint32_t max_size,
                   uint32_t size,
                   uint32_t min,
                   uint32_t max,
                   s_align* p) : max_size(max_size), size(size), min_index(min), max_index(max), ptr(p) {}
};

/*! @fn compute_lis_alignment()
    @brief compute the LIS and sequence alignment for a read
    @param uint32_t size_ambiguous_nt
    @param uint32_t readhit
    @param vector< id_win >& id_win_hits
    @param kmer_origin* positions_tbl
    @param uint16_t* read_max_SW_score
    @param bool& search
    @param int32_t* best_x
    @param uint64_t readn
    @param int32_t* num_alignments_x
    @param uint32_t readlen
    @param uint32_t lnwin_index_num
    @param uint16_t index_num
    @param uint64_t* reference_seq_len
    @param char* myread
    @param int32_t* ambiguous_nt
    @param int8_t* mat
    @param char** reference_seq
    @param long gap_open
    @param long gap_extension
    @param uint32_t minimal_score_index_num
    @param vector<bool>& read_hits
    @param uint64_t& total_reads_mapped
    @param vector<uint64_t>& reads_matched_per_db
    @param uint16_t part
    @param map<uint64_t, alignment_struct >& read_hits_align_info
    @param uint32_t max_SW_score
    @param bool& read_to_count
    @param uint64_t& total_reads_mapped_cov
    @param vector<bool>& read_hits_denovo
    @param char filesig
    @param uint64_t strs
    @param uint32_t file_s
    @param uint32_t file_sections
    @param char** reads
    @param char* finalnt
    @param double gumbel_lambda_index_num
    @param double gumbel_K_index_num
    @param uint64_t full_ref_index_num
    @param uint64_t full_read_index_num
    @param ofstream& acceptedblast
    @param ofstream& acceptedsam
    @return none
*/
void
compute_lis_alignment(uint32_t size_ambiguous_nt /**< number of ambiguous nucleotides in a read */,
                      uint32_t readhit /**< number of seeds matches between read and database */,
                      vector< id_win >& id_win_hits /**< mini burst trie seed IDs for current window on read */,
                      kmer_origin* positions_tbl /**< mini burst trie seed IDs pointing to (L+1)-mer positions in reference database */,
                      uint16_t* read_max_SW_score /**< array storing number of times a read was aligned with maximum alignment score */,
                      bool& search /**< if true, continue searching for seeds at lower granularity intervals across read */,
                      int32_t* best_x /**< array storing number of candidate reference sequences to explore for read (based on min_lis) */,
                      uint64_t readn /**< read number */,
                      int32_t* num_alignments_x /**< array storing number of alignments output for read */,
                      uint32_t readlen /**< length of read */,
                      uint32_t lnwin_index_num /**< seed length (L) used to index reference database */,
                      uint16_t index_num /**< reference database number */,
                      uint64_t* reference_seq_len /**< array storing lengths of all reference sequences in an index */,
                      char* myread /**< read sequence */,
                      int32_t* ambiguous_nt /**< array storing positions of ambiguous nucleotides in read */,
                      int8_t* mat /**< Smith-Waterman scoring matrix */,
                      char** reference_seq /**< array storing raw nucleotide reference sequences in an index */,
                      long gap_open /**< Smith-Waterman gap open score */,
                      long gap_extension /**< Smith-Waterman gap extension score */,
                      uint32_t minimal_score_index_num /**< minimal Smith-Waterman score allowed based on E-value */,
                      vector<bool>& read_hits /**< array storing a boolean for each read, 1 if a read was aligned, 0 otherwise */,
                      uint64_t& total_reads_mapped /**< total number of reads aligned passing E-value threshold */,
                      vector<uint64_t>& reads_matched_per_db /**< number of reads aligned per reference database (ex. 16S bacteria and 18S eukarya) */,
                      uint16_t part /**< index part number */,
                      map<uint64_t, alignment_struct >& read_hits_align_info /**< a map storing alignment information for each aligned read */,
                      uint32_t max_SW_score /**< maximum achievable alignment score for read */,
                      bool& read_to_count /**< flag to count only one alignment per read */,
                      uint64_t& total_reads_mapped_cov /**< total number of reads mapped passing E-value threshold & %id and/or %query coverage thresholds */,
                      vector<bool>& read_hits_denovo /**< array toring a boolean for each read, 1 if > %id and > %coverage threshold, 0 otherwise */,
                      char filesig /**< @ or > depending on FASTQ or FASTA reads input file */,
                      uint64_t strs /**< number of reads in input file */,
                      uint32_t file_s /**< file section currently memory mapped */,
                      uint32_t file_sections /**< total number of file sections to be memory mapped */,
                      char** reads /**< array storing pointers to sequences in the memory mapped buffer */,
                      char* finalnt /**< final character in memory map */,
                      double gumbel_lambda_index_num /**< Lambda value for current index */,
                      double gumbel_K_index_num /**< K value for current index */,
                      uint64_t full_ref_index_num /**< size of all reference sequences in current index */,
                      uint64_t full_read_index_num /**< size of all read sequences in current memory map */,
                      ofstream& acceptedblast /**< output stream for aligned reads in BLAST format */,
                      ofstream& acceptedsam /**< output stream for aligned reads in SAM format */);

#endif