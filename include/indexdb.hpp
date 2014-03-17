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
 * file: indexdb.hpp
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */


#ifndef INDEXDB_H
#define INDEXDB_H

#include <unistd.h>  /// lseek
#include <fcntl.h>   /// file checking, open/close
#include <cstdlib>
#include <cstdio>
#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <sys/mman.h> 
#include <set>
#include <sys/time.h>
#include <math.h>
#include <deque>

using namespace std;

extern timeval t;

#define TIME(x) gettimeofday(&t, NULL); x = t.tv_sec + (t.tv_usec/1000000.0);

/// maximum number of bytes in a bucket prior to bursting, should be <= L1 cache size and a power of 2 for efficiency
   #define THRESHOLD 128

extern bool verbose;

/// output for verbose option
#define eprintf(format, ...) do { if (verbose) fprintf(stderr, format, ##__VA_ARGS__);} while(0)


/* use ascii decimal value of a letter as an offset in this array such that:
		A/a -> 0, C/c -> 1, G/g -> 2, T/t -> 3, U/u -> 3
		R/r -> 0, Y/y -> 1, S/s -> 2, W/w -> 1, K/k -> 2
		M/m -> 0, B/b -> 1, D/d -> 0, H/h -> 0, V/v -> 0
		N/n -> 0	 */

extern const char map_nt[122];


struct NodeElement
{
	/// a pointer to a bucket or another trie node
	union 
	{
	   void* bucket;
	   NodeElement* trie;
	} whichnode;

	/// current size (in bytes) of bucket
	unsigned int size; 

	/// flags whether the node leads to a child trie node, is empty or leads to a bucket node;
	///    0 :: empty;
	///    1 :: trie node;
 	///    2 :: bucket
	char flag;
};

/// the sequence number and index position at which a 19-mer exists on the reference sequence; these values *must* be positive
struct seq_pos
{
	uint32_t pos; /// position on the sequence
	uint32_t seq; /// the sequence
};

struct kmer_origin
{
	seq_pos* arr; /// pointer to seq_pos array
	uint32_t size; /// number of 19-mer occurrences
};

struct kmer
{
	NodeElement* trie_F; /// pointer to forward mini burst trie
	NodeElement* trie_R; /// pointer to reverse mini burst trie
	uint32_t count; /// count of 9-mers
};

/// data structure to store information on index parts built
struct index_parts_stats {
    unsigned long int start_part; /// where the section starts in the file
    unsigned long int seq_part_size; /// number of bytes of reference sequences to read
    uint32_t numseq_part; /// the number of sequences in this part 
};



#endif
