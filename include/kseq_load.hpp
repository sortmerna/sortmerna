/*
 * @file kseq_load.hpp
 * @brief Function and variable definitions for kseq_load.cpp
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
 * along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
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
 
#ifndef KSEQ_LOAD_H
#define KSEQ_LOAD_H

#include "common.hpp"
#include "kseq.h"
#ifdef HAVE_LIBZ
#include <zlib.h>
KSEQ_INIT(gzFile, gzread)
#else
KSEQ_INIT(int, read)
#endif



/*! @fn load_reads()
    @brief load reads into buffer using kseq library
    @param char* inputreads
    @param char* raw
    @param uint64_t number_total_read
    @param off_t full_file_size
    @param char*& finalnt
    @return char** reads
    @version Feb 08, 2016 
*/
char**
load_reads(char* inputreads,
           char*& raw,
           uint64_t number_total_read,
           off_t full_file_size,
           char*& finalnt);

#endif