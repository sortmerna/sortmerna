/**
 * @file paralleltraversal.hpp
 * @brief Function and variable definitions for paralleltraversal.cpp
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

#ifndef PARALLELTRAVERSAL_H
#define PARALLELTRAVERSAL_H

#include <iomanip>
 
#include "outputformats.hpp"
#include "load_index.hpp"
#include "traverse_bursttrie.hpp"
#include "alignment.hpp"
#include "mmap.hpp"
#include "kseq_load.hpp"


/*! @fn check_file_format()
    @brief check input reads file format (FASTA, FASTQ or unrecognized)
    @param char* inputreads
    @param char& filesig
    @return bool
    @version Feb 15, 2016
 */
bool
check_file_format(char* inputreads /**< pointer to query reads file */,
                  char& filesig /**< first character of sequence label */);

/*! @fn compute_read_stats()
    @brief compute total number of reads in file and their combined length
    @param char* inputreads
    @param uint64_t& number_total_read
    @param uint64_t& full_read_main
    @param off_t& full_file_size
    @version Feb 15, 2016
 */
void compute_read_stats(char* inputreads /**< pointer to query reads file */,
                        uint64_t& number_total_read /**< total number of reads */,
                        uint64_t& full_read_main /**< total number of nucleotides in all reads */,
                        off_t& full_file_size /**< the size of the full reads file (in bytes) */);


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

    @param char* inputreads
    @param bool have_reads_gz
    @param *ptr_filetype_ar 
    @param *ptr_filetype_or
    @param long match
    @param long mismatch
    @param long gap_open
    @param long gap_extension
    @param long score_N
    @param vector< vector<uint32_t> > 
    @param int argc
    @param char **argv
    @param bool yes_SQ
    @param vector< pair<string,string> >& myfiles
    @return void
    @version Feb 10, 2016 
*/
void
paralleltraversal (char* inputreads /**< pointer to query reads file */,
                   bool have_reads_gz /**< if true, input reads file is in compressed format */,
                   char* ptr_filetype_ar /**< pointer to string for aligned seqeunces filepath */,
                   char* ptr_filetype_or /**< pointer to string for rejected sequences filepath */,
                   long match /**< SW match reward score (positive) */,
                   long mismatch /**< SW mismatch penalty score (negative) */,
                   long gap_open /**< SW gap open penalty score (positive) */,
                   long gap_extension /**< SW gap extend penalty score (positive) */,
                   long score_N /**< SW penalty for ambiguous nucleotide (negative) */,
                   vector< vector<uint32_t> >& skiplengths /**< skiplengths, three intervals at which to place seeds on read */,
                   int argc /**< number of arguments passed to SortMeRNA */,
                   char **argv /**< argument string passed to SortMeRNA */,
                   bool yes_SQ /**< if true, include @SQ tags in SAM output */,
                   vector< pair<string,string> >& myfiles /**< vector of (FASTA file, index name) pairs for loading index */,
                   bool exit_early /**< if true, exit program if reads file is not FASTA or FASTQ, or reads files or reference file is empty */);

#endif //~parallel_traversal_h
