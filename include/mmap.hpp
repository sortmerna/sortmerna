/**
 * @file mmap.hpp
 * @brief Load data (reads) using mmap.
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

#ifndef MMAP_H
#define MMAP_H

#include "load_index.hpp"

/*! @fn mmap_reads()
    @brief mmap file
    @param off_t partial_file_size
    @param char* inputreads
    @param off_t offset_map
    @param char* raw
    @param char filesig
    @param uint32_t file_s
    @param uint32_t file_sections
    @param int32_t offset_pair_from_top
    @param char* split_read_ptr
    @param char* split_read
    @param uint64_t &strs
    @param char*& finalnt
    @param uint32_t &reads_offset_f
    @param uint32_t &reads_offset_e
    @return char**
    @version Feb 08, 2016 
*/
char**
mmap_reads(off_t partial_file_size /**< number of bytes in memory map buffer */,
           char* inputreads /**< pointer to query reads file */,
           off_t offset_map /**< the offset from the start of the reads file for mmap */,
           char*& raw /**< the address of the new mapping from mmap */,
           char filesig /**< @ or > to identify FASTQ or FASTA file */,
           uint32_t file_s /**< current file section number */,
           uint32_t file_sections /**< number of file sections (size of reads file divided by mmap buffer size) */,
           int32_t &offset_pair_from_top /**< number of lines to offset at the top of the current file section */,
           char* split_read_ptr /**< pointer to the position in the split read where to attach the connecting part of the split read (and possibly its pair) */,
           char* split_read /**< pointer to the split read (the read which is split between any two file sections) */,
           uint64_t &strs /**< number of reads in current memory mapped section */,
           char*& finalnt /**< pointer to final character in memory mapped buffer */,
           uint32_t &reads_offset_f /**< the length of the split read in file part i+1 (from beginning of file part) */,
           uint32_t &reads_offset_e /**< the length of the split read in file part i (from end of file part) */,
           uint32_t min_lnwin /**< the minimum seed length used to index reference databases */);


/*! @fn unmmap_reads()
    @brief unmap file
    @param char*& raw
    @param off_t partial_file_size
    @return void
    @version Feb 08, 2016 
*/
void
unmmap_reads(char*& raw /**< the address of the new mapping from mmap */,
             off_t partial_file_size /**< number of bytes in memory map buffer */);

#endif
