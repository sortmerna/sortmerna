/**
 * @file load_index.hpp
 * @brief header file for load_index.cpp
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

#ifndef LOADINDEX_H
#define LOADINDEX_H

#include <string>
#include <vector>

#include "../alp/sls_alignment_evaluer.hpp" //! ALP program for computing the Gumbel parameters

#include "indexdb.hpp"

using namespace Sls;


/*! @brief Map nucleotides to integers.
Ambiguous letters map to 4.
{A/a,C/c,G/g,T/t,U/u} = {0,1,2,3,3} respectively.
*/
const char nt_table[128] = {
	'4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',
	'4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',
	'4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',
	'4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',
	'4', '0', '4', '1',  '4', '4', '4', '2',  '4', '4', '4', '4',  '4', '4', '4', '4',
	'4', '4', '4', '4',  '3', '3', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4',
	'4', '0', '4', '1',  '4', '4', '4', '2',  '4', '4', '4', '4',  '4', '4', '4', '4',
	'4', '4', '4', '4',  '3', '3', '4', '4',  '4', '4', '4', '4',  '4', '4', '4', '4'
};


 /*! @fn load_index_stats()
    @brief Load reference database index statistics.
    @detail
    @param vector< pair<string,string> >& myfiles
    @param char** argv
    @param int argc
    @param bool yes_SQ
    @param char* acceptedstrings_sam
    @param long _match
    @param long _mismatch
    @param long _gap_open
    @param long _gap_extension
    @param vector<vector<uint32_t> >&
    @param vector<uint16_t>& num_index_parts
    @param vector<vector<index_parts_stats> >& index_parts_stats_vec
    @param vector<uint64_t>& full_ref
    @param vector<uint64_t>& full_read
    @param vector<uint32_t>& lnwin
    @param vector<uint32_t>& partialwin
    @param vector<uint32_t>& minimal_score 
    @param uint64_t number_total_read
    @param vector<pair<double, double> >& gumbel
    @param vector<uint64_t>& numbvs 
    @param vector<uint64_t>& numseq
    @return none
*/
void
load_index_stats(
	std::vector<std::pair<std::string, std::string>> & myfiles, /**< vector of (FASTA file, index name) pairs for loading index */
    char** argv, /**< command line for executing SortMeRNA */
    int argc, /**< number of arguments in command line for executing SortMeRNA */
    bool yes_SQ, /**< if true, include @SQ tags in SAM output */
    char* acceptedstrings_sam, /**< pointer to output SAM file */
    long _match, /**< Smith-Waterman score for a match */
    long _mismatch, /**< Smith-Waterman score for a mismatch */
    long _gap_open, /**< Smith-Waterman score for gap opening */
    long _gap_extension, /**< Smith-Waterman score for gap extension */
    std::vector<std::vector<uint32_t>> & skiplengths, /**< skiplengths, three intervals at which to place seeds on read */
    std::vector<uint16_t>& num_index_parts, /**< number of index files */
    std::vector<std::vector<index_parts_stats>> & index_parts_stats_vec, /**< statistics for index files */
    std::vector<uint64_t>& full_ref, /**< corrected size of each reference index (for computing E-value) */
    std::vector<uint64_t>& full_read, /**< corrected size of reads (for computing E-value) */
    std::vector<uint32_t>& lnwin, /**< length of seed (sliding window L) */
    std::vector<uint32_t>& partialwin, /**< length of seed/2 */
    std::vector<uint32_t>& minimal_score, /**< minimal SW score in order to reach threshold E-value */
    uint64_t number_total_read, /**< total number of reads in input reads file */
    std::vector<pair<double, double> >& gumbel, /**< Gumbel parameters Lambda and K */
    std::vector<uint64_t>& numbvs, /**< number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward */
    std::vector<uint64_t>& numseq /**< total number of reference sequences in one complete reference database */
);

 /*! @fn load_index()
    @brief Load reference database index.
    @detail
    @param char* ptr_dbindex
    @param string part_str
    @param kmer*& lookup_tbl
    @param kmer_origin*& positions_tbl
    @param uint32_t& number_elements
    @param uint32_t lnwin
    @return none
*/
void
load_index(
	char* ptr_dbindex, /**< pointer to index file name */
    std::string part_str, /**< index part number */
    kmer*& lookup_tbl, /**< reference to L/2-mer look up table */
    kmer_origin*& positions_tbl, /**< reference to (L+1)-mer positions table */
    uint32_t& number_elements, /**< number of positions in (L+1)-mer positions table */
    uint32_t lnwin /**< length of seed (sliding window L) */
);

 /*! @fn load_ref()
    @brief Load reference database.
    @detail
    @param char* ptr_dbfile
    @param char* buffer
    @param char** reference_seq
    @param uint64_t* reference_seq_len
    @param uint64_t seq_part_size
    @param uint64_t numseq_part
    @param uint64_t start_part
    @param bool load_for_search
    @return none
*/
void
load_ref(
	char* ptr_dbfile, /**< pointer to reference database file */
    char* buffer, /**< pointer to memory slot for storing reference database */
    char** reference_seq, /**< array of pointers to sequences in buffer */
    uint64_t* reference_seq_len, /**< array of lengths for each sequence in buffer */
    uint64_t seq_part_size, /**< size of memory to allocate for buffer */
    uint64_t numseq_part, /**< number of sequences in part of database indexed */
    uint64_t start_part, /**< index of first sequence in current index */
    bool load_for_search /**< if true, compute sequence length; if false, only load sequence */
);

#endif