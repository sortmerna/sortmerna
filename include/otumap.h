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

#pragma once

#include <string>
#include <vector>
#include <map>

// forward
struct Runopts;
class Refstats;
struct Readstats;
class Readfeed;
class References;
class KeyValueDatabase;

class OtuMap {
	// Clustering of reads around references by similarity i.e. {ref: [read, read, ...] , ref : [read, read...] , ...}
	// calculated after alignment is done on all reads
	//  TODO: Store in DB ? Can be very big.
	std::vector<std::map<std::string, std::vector<std::string>>> mapv;
	//std::map<std::string, std::vector<std::string>> otu_map;
	//          |_Ref_ID       |_Read_IDs
public:
	std::filesystem::path fmap;
	uint64_t total_otu; // total count of OTU groups in otu_map
public:
	OtuMap(int numThreads=1);
	void push(int idx, std::string& ref_seq_str, std::string& read_seq_str);
	void merge();
	void write();
	void init(Runopts& opts);
	size_t count_otu();
};

void fill_otu_map(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts, bool is_write=true);
