/*
 * @file triestats.hpp
 * @brief Compute burst trie statistics.
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2012-17 Bonsai Bioinformatics Research Group
 * @copyright 2014-17 Knight Lab, Department of Pediatrics, UCSD, La Jolla
 * @copyright 2016- Clarity Genomics Inc
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

#ifndef triestats_h
#define triestats_h

#include "common.hpp"
#include "bursttrie.hpp"
#include <fstream>

using namespace std;

extern int total_num_trie_nodes;
extern size_t size_of_all_buckets;
extern int total_num_buckets;

void traverse_trie( NodeElement* root, int depth);

#endif
