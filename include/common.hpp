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
 * file: common.hpp
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */

#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <bitset>
#include <vector>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sys/time.h>
#include <set>
#include <stdint.h>

/// length of reads
#define READLEN 30000

/// detect chimeric reads (under development)
//#define chimera


/// Debugging memory map and splitting of reads file into multiple file
/// sections (output piped to stdout)

//#define debug_mmap

/// Debugging the collection of seeds for constructing an LIS and the
/// Smith-Waterman alignment using the LIS (output piped to stdout)

//#define debug_align

/// Debugging for seeing the output for writing mini-burst tries into binary
/// during indexdb_rna and also for reading from binary to re-construct the
/// mini-burst tries in sortmerna (both will be actived with this define statement)
/// Outputs for both must be equal (piped to stdout)
/// Will require C++11 compiler suppport, add -std=c++0x to CXXFLAGS & CFLAGS in Makefile

//#define see_binary_output

/// minimum number of hits of a read on an rrna database required for output
extern int32_t numcpu_gv;
extern bool verbose;

/// flag to output accepted reads by matching database
extern bool bydbs_gv;

extern bool pairedin_gv;

extern bool pairedout_gv;

/// output chimeric reads
extern bool chimeraout_gv;

/// output overall statistics file
extern bool logout_gv;

/// output FASTA/FASTQ reads passing E-value threshold but having < %id and < %coverage scores for de novo OTU construction
extern bool de_novo_otu_gv;

/// forward and reverse strands to search
extern bool forward_gv;
extern bool reverse_gv;

/// output OTU map
extern bool otumapout_gv;

/// flag to include pid in output file names
extern bool pid_gv;

/// size of partial section of reads file to mmap
extern long unsigned int map_size_gv;

/// SAM output
extern bool samout_gv;

/// BLAST-like output
extern bool blastout_gv;

/// BLAST-like output format
extern int32_t blast_outfmt;

/// FASTA/Q output
extern bool fastxout_gv;

/// output first, best or all alignments
extern int32_t min_lis_gv;

/// number of best alignments per read to output
extern int16_t num_best_hits_gv;

/// output first num_alignments_gv alignments
extern int32_t num_alignments_gv;

/// number of seed hits before searching for candidate LCS
extern int32_t seed_hits_gv;

/// number of nucleotides to add to each edge of the alignment region prior to extension
extern int32_t edges_gv;

/// flag to turn off heuristic for stopping index search after finding a 0-error match, instead collect all 0-error and 1-error matches
extern bool full_search_gv;

/// E-value threshold
extern double evalue;

/// minimum %id to keep alignment
extern double align_id;

/// minimum %coverage to keep alignment
extern double align_cov;

/// interpret option --edges as percent
extern bool as_percent_gv;

/// do not mask low-occurring (L/2)-mers for seeds of length L
extern bool nomask_gv;

/// if true, print all reads in SAM and/or BLAST output (aligned and non-aligned)
extern bool print_all_reads_gv;

/// output for verbose option
#define eprintf(format, ...) do { if (verbose) fprintf(stderr, format, ##__VA_ARGS__);} while(0)

#endif

