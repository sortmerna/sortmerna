/**
 * @file load_index.hpp
 * @brief header file for traverse_bursttrie.cpp
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2013-16 Bonsai Bioinformatics Research Group
 * 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
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

#ifndef LOADINDEX_H
#define LOADINDEX_H

#include <deque>
#include "indexdb.hpp"
 //! ALP program for computing the Gumbel parameters
#include "../alp/sls_alp_data.hpp"
#include "../alp/sls_alp_sim.hpp"
#include "../alp/gumbel_params.hpp"

using namespace std;

extern char nt_table[128];

 /*! @fn load_index_stats()
    @brief Load reference database index statistics.
    @detail
    @param vector< pair<string,string> >& myfiles, vector of (FASTA file, index name) pairs for loading index 
    @param char** argv, command line for executing SortMeRNA
    @param int argc, number of arguments in command line for executing SortMeRNA
    @param bool yes_SQ, boolean for outputting SQ fields in SAM output
    @param char* acceptedstrings_sam, pointer to output SAM file
    @param int32_t _match, Smith-Waterman score for a match
    @param int32_t _mismatch, Smith-Waterman score for a mismatch
    @param int32_t _gap_open, Smith-Waterman score for gap opening
    @param int32_t _gap_extension, Smith-Waterman score for gap extension
    @param vector<vector<uint32_t> >& skiplengths, three intervals at which to place seeds on read
    @param vector<uint16_t>& num_index_parts, number of index files
    @param vector<vector<index_parts_stats> >& index_parts_stats_vec, statistics for index files
    @param vector<uint64_t>& full_ref, corrected size of each reference index (for computing E-value)
    @param vector<uint64_t>& full_read, corrected size of reads (for computing E-value)
    @param vector<uint32_t>& lnwin, length of seed (sliding window)
    @param vector<uint32_t>& partialwin, length of seed/2
    @param vector<uint32_t>& minimal_score, minimal SW score in order to reach threshold E-value
    @param uint32_t number_total_read, total number of reads in input reads file
    @param vector<pair<double, double> >& gumbel, Gumbel parameters Lambda and K
    @param vector<uint32_t>& numbvs, number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward
    @param vector<uint32_t>& numseq, total number of reference sequences in one complete reference database
    @return none
*/
void
load_index_stats(vector< pair<string,string> >& myfiles,
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
                 vector<uint32_t>& numseq );

 /*! @fn load_index()
    @brief Load reference database index.
    @detail
    @param char* ptr_dbindex, 
    @param string part_str,
    @param kmer*& lookup_tbl,
    @param kmer_origin*& positions_tbl,
    @param uint32_t& number_elements,
    @param uint32_t lnwin
    @return none
*/
void
load_index(char* ptr_dbindex,
           string part_str,
           kmer*& lookup_tbl,
           kmer_origin*& positions_tbl,
           uint32_t& number_elements,
           uint32_t lnwin);

 /*! @fn load_ref()
    @brief Load reference database.
    @detail
    @param char* ptr_dbfile, 
    @param char* buffer,
    @param char** reference_seq,
    @param uint32_t* reference_seq_len,
    @param uint64_t seq_part_size,
    @param uint32_t numseq_part,
    @param uint64_t start_part,
    @param bool load_for_search,
    @return none
*/
void
load_ref(char* ptr_dbfile,
         char* buffer,
         char** reference_seq,
         uint32_t* reference_seq_len,
         uint64_t seq_part_size,
         uint32_t numseq_part,
         uint64_t start_part,
         bool load_for_search);
#endif