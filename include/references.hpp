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
 * FILE: references.hpp
 * Created: Nov 06, 2017 Mon
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <algorithm>

#include "common.hpp" // Format, FASTA_HEADER_START, FASTQ_HEADER_START

// forward
class Refstats;
struct Runopts;

class References {
public:
	struct BaseRecord
	{
		size_t nid; // position of the sequence in the Reference file [0...'number of sequences in the ref.file - 1']
		std::string id; // ID from header
		std::string header;
		std::string sequence;
		std::string quality; // "" (fasta) | "xxx..." (fastq)
		BIO_FORMAT format; // FASTA | FATSQ
		bool isEmpty;
		BaseRecord(): isEmpty(true) {}
		void clear()
		{
			header.clear();
			sequence.clear();
			quality.clear();
			isEmpty = true;
		}

		std::string getId() {
			std::string id = header.substr(0, header.find(' '));
			id.erase(id.begin(), std::find_if(id.begin(), id.end(), [](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));
			return id;
		}
	};

	std::vector<BaseRecord> buffer; // Container for references TODO: change name?

	References(): num(0), part(0) {}
	//~References() {}

	void load(uint32_t idx_num, uint32_t idx_part, Runopts & opts, Refstats & refstats); // load references into the buffer given index number and index part
	void convert_fix(std::string & seq); // convert sequence to numberical form and fix ambiguous chars
	std::string convertChar(int idx); // convert numerical form to char string
	/*
	* For debugging needs.
	* find a reference index given a header.
	*/
	int findref(const std::string& id);
	void unload();

public:
	uint16_t num; // number of the reference file currently loaded
	uint16_t part; // part of the reference file currently loaded

//private:
//	bool load_for_search;
}; // ~class References
