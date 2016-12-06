/*
 * @file indexdb.hpp
 * @brief Function and variable definitions for indexdb.cpp
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2012-16 Bonsai Bioinformatics Research Group
 * @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
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
 * along with SortMeRNA.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 *
 */

 
#ifndef INDEXDB_H
#define INDEXDB_H

#include <sys/types.h>
#include "ssw.h"
#include "common.hpp"


/*! @brief The size of an entry in a burst trie bucket 

    Each bucket entry is of ENTRYSIZE bytes. The first
    4 bytes are to store the seed (ex. 18-mer) and the
    last 4 bytes are to store the index for the 
    positions lookup table.<br/>

    The seed is encoded as 2 bits/nt, such that
    16 nucleotides can fit into a 4 byte integer.
*/ 
#define ENTRYSIZE (2*sizeof(uint32_t))

// maximum number of bytes in a bucket prior to bursting, should be <= L1 cache size and a power of 2 for efficiency
#define THRESHOLD 128

struct NodeElement
{
	// a pointer to a bucket or another trie node
	union 
	{
	   void* bucket;
	   NodeElement* trie;
	} whichnode;

	// current size (in bytes) of bucket
	unsigned int size; 

	// flags whether the node leads to a child trie node, is empty or leads to a bucket node;
	//    0 :: empty;
	//    1 :: trie node;
 	//    2 :: bucket
	char flag;
};

// the sequence number and index position at which a 19-mer exists on the reference sequence; these values *must* be positive
struct seq_pos
{
	uint32_t pos; // position on the sequence
	uint32_t seq; // the sequence
};

struct kmer_origin
{
	seq_pos* arr; // pointer to seq_pos array
	uint32_t size; // number of 19-mer occurrences
};

struct kmer
{
	NodeElement* trie_F; // pointer to forward mini burst trie
	NodeElement* trie_R; // pointer to reverse mini burst trie
	uint32_t count; // count of 9-mers
};

// data structure to store information on index parts built
struct index_parts_stats {
    unsigned long int start_part; // where the section starts in the file
    unsigned long int seq_part_size; // number of bytes of reference sequences to read
    uint32_t numseq_part; // the number of sequences in this part
};



#endif
