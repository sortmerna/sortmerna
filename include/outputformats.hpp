/*
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * Copyright (C) 2012-2014 Bonsai Bioinformatics Research Group
 *
 * OTU-picking extensions developed in the Knight Lab, BioFrontiers Institute,
 * University of Colorado at Boulder, Boulder, CO
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
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */

/** @file outputformats.hpp */ 
 
#ifndef OUTPUTFORMATS_H
#define OUTPUTFORMATS_H

#include <istream>

//! SIMD Smith-Waterman alignment library
#include "indexdb.hpp"

using namespace std;

void
report_blast (ofstream &fileout,
                   s_align* a,
                   char* read_name,
                   char* read_seq,
                   char* read_qual,
                   char* ref_name,
                   char* ref_seq,
                   double evalue,
                   uint32_t readlen,
                   uint32_t bitscore,
                   bool strand,
                   double id,
                   double coverage,
                   uint32_t mismatches,
                   uint32_t gaps);


void
report_sam (ofstream &fileout,
                 s_align* a,
                 char* read_name,
                 char* read_seq,
                 char* read_qual,
                 char* ref_name,
                 char* ref_seq,
                 uint32_t readlen,
                 bool strand,
                 uint32_t diff);

void
report_fasta (char* acceptedstrings,
                   char* ptr_filetype_or,
                   char* ptr_filetype_ar,
                   char** reads,
                   int32_t strs,
                   vector<bool>& read_hits,
                   uint32_t file_s,
                   char* finalnt);

void
report_denovo(char *denovo_otus_file,
                   char **reads,
                   int32_t strs,
                   vector<bool>& read_hits_denovo,
                   uint32_t file_s,
                   char *finalnt);

void report_biom (char* biomfile);

#endif