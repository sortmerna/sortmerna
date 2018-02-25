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
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 *
 */

 /** @file common.hpp */

#ifndef COMMON_H
#define COMMON_H

#include <sys/time.h>
#include "config.h"

const char FASTA_HEADER_START = '>';
const char FASTQ_HEADER_START = '@';
//const char WIN_ENDLINE = '\r';
//const char UNI_ENDLINE = '\n';
const int QUEUE_SIZE_MAX = 10; // TODO: move to process options?

enum class Format { FASTA, FASTQ }; // format of Reads and References files. Used in References and Read
enum class BlastFormat { TABULAR, REGULAR}; // format of the Blast output
enum class CompressionType { XPRESS, ZLIB };

/*! @brief Map nucleotides to integers.
Ambiguous letters map to 4.
{A/a,C/c,G/g,T/t,U/u} = {0,1,2,3,3} respectively.
*/
const char nt_table[128] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

extern timeval t;

/*! @brief Macro for timing */
#define TIME(x) gettimeofday(&t, NULL); x = t.tv_sec + (t.tv_usec/1000000.0);

/*! @brief Print function for verbose mode */
#define eprintf(format, ...) do {if (false) fprintf(stdout, format, ##__VA_ARGS__);} while(0)

/*! @brief start color text red */
#if defined(_WIN32)
#define RED    ""
#define GREEN  ""
#define YELLOW ""
#define BOLD   ""
#define UNDL   ""
#define COLOFF ""
#else
#define RED    "\033[0;31m"
#define GREEN  "\033[0;32m"
#define YELLOW "\033[0;33m\"
#define BOLD   "\033[1m"
#define UNDL   "\033[4m" // underline
#define COLOFF "\033[0m" // color off
#endif



/*! @brief Maximum length of input reads
	(not limited to this length algorithmically)
*/
#define READLEN 30000

#endif

