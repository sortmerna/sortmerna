#pragma once
/**
* FILE: references.hpp
* Created: Nov 06, 2017 Mon
*/

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
		std::string header;
		std::string sequence;
		std::string quality; // "" (fasta) | "xxx..." (fastq)
		Format format; // FASTA | FATSQ
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
	~References() {}

	void load(uint32_t idx_num, uint32_t idx_part, Runopts & opts, Refstats & refstats); // load references into the buffer given index number and index part
	void convert_fix(std::string & seq); // convert sequence to numberical form and fix ambiguous chars
	std::string convertChar(int idx); // convert numerical form to char string
	void clear();

public:
	uint16_t num; // number of the reference file currently loaded
	uint16_t part; // part of the reference file currently loaded

private:
	bool load_for_search;
}; // ~class References
