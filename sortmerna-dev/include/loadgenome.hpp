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
 * file: loadgenome.hpp
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */

#include "paralleltraversal.hpp"
#include <queue>

#ifndef minoccur_h
#define minoccur_h

/*
 *
 * FUNCTION 	: void find_minoccur   ( int minlenread )
 * PURPOSE	: find the minimum possible occurrence of a k-mer in the database for which a parallel
 *		  traversal of the universal Levenshtein automaton and the Burst trie is carried out 	  
 *
 **************************************************************************************************************/
void find_minoccur ( int );

#ifdef WINDOWS
/*
 *
 * FUNCTION 	: int isnewline   ( char *ch )
 * PURPOSE	: parse '\n', '\r' or '\r\n' newlines
 * PARAMETERS	: pointer to character being observed	  
 *
 **************************************************************************************************************/
int isnewline( char *ch );
#endif

/* rRNA database file (fasta) */
extern char *rrnadbfile;

/* the minimum number of occurrences of an s/2-mer hash id in the rrna database */
extern unsigned int MINOCCUR;

/* maximum/minimum length of an rRNA sequence in the database */
extern int maxlen;
extern int minlen;

/* the number of rRNAs in the rRNA database */
extern int rrnastrings;

extern unsigned int mask32;
extern char map_nt[122];

extern hashid *kmerf;
extern hashid *kmerr;


#endif

