#pragma once
/**
 * @file paralleltraversal.hpp
 * @brief Function and variable definitions for paralleltraversal.cpp
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

#include "stdint.h"
#include <string>
#include <vector>
 
#include "options.hpp"

// forward
struct Readstats;
class Output;

/*! @fn align()
	@brief Traverse the query input and indexed database and output
		   alignments passing the E-value threshold
	@detail The following methods will be executed:
	<ol>
	  <li> compute the gumbel parameters (lamda and K) using ALP,
		   load the index fully or in parts (depending on how
		   it was built) </li>
	  <li> using 3 intervals, scan over the read and collect all
		   L-mers on the read which match to the reference
		   index with at most 1 error. This is done using
		   parallel traversal between the index and the
		   Levenshtein automaton </li>
	  <li> if enough L-mers were collected, extend them into
		   longer matches using the Longest Increasing
		   subsequence (LIS) of positions where the L-mers
		   matched on the read and the associated reference
		   sequences </li>
	  <li> if the LIS is long enough, use the starting positions
		   of the LIS to estimate the starting position
		   of an alignment and pass this reference segment and
		   read to SSW </li>
	  <li> if the alignment score is at least the minimum score
		   corresponding to the E-value threshold, keep the read,
		   otherwise continue searching for other LIS or more
		   L-mers using smaller intervals </li>
	</ol>
*/
void align(Runopts & opts, Readstats & readstats, Output & output);

// ~PARALLELTRAVERSAL_H
