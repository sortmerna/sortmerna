/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

/*
 * file: summary.hpp
 *
 * post-alignment tasks like calculating statistics 
 */

#pragma once

#include <string>
#include <vector>

// forward
struct Runopts;
class Refstats;
struct Readstats;
class Readfeed;
class References;
class KeyValueDatabase;

void writeSummary(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts);

/**
 * Summary report (log) data structure
 */
class Summary {
public:
	// vars
	bool is_de_novo;
	bool is_otumapout;
	std::string cmd;
	std::string timestamp;
	std::string pid_str;
	uint64_t total_reads;
	uint64_t total_mapped; // total passing SW
	uint64_t total_denovo; // total SW + !ID + !COV
	uint64_t total_id_cov; // total SW + ID + COV
	uint64_t total_otu; // total number of OTU groups
	uint32_t min_read_len;
	uint32_t max_read_len;
	uint64_t all_reads_len;
	std::vector<std::pair<std::string, float>> db_matches;

	// methods
	Summary();
	void write(Refstats& refstats, Readstats& readstats, Runopts& opts);
	std::string to_string(Refstats& refstats, Runopts& opts);
};
