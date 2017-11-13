#pragma once
/**
* FILE: references.hpp
* Created: Nov 06, 2017 Mon
*/

#include <vector>

#include "options.hpp"
#include "index.hpp"
#include "common.hpp"

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
	};

	References(Runopts & opts, Index & index) : opts(opts), index(index) {}
	~References() {}

	void load(uint32_t idx_num, uint32_t idx_part); // load refrences into the buffer given index number and index part
	void convert_fix(std::string & seq); // convert sequence to numberical form and fix ambiguous chars

	std::vector<BaseRecord> buffer; // Container for references TODO: change name?
private:
	Runopts & opts;
	Index & index;
	//uint16_t idx_num; // number of currently loaded index/reference file - Index::index_num
	//	char* ptr_dbfile; // opts
	// use vector instead of 3 variables below
	//char* buffer; // holds sequences from a part of a reference file. Calculated after Index::load_stats is called.
	//char** reference_seq; // pointers to the start of every sequence in the '*buffer'
	//uint64_t* reference_seq_len; // lengths of the sequences in the '*buffer'

	//uint64_t seq_part_size; // take directly from Index::index_parts_stats_vec
	// uint64_t numseq_part; // take directly from Index::index_parts_stats_vec
	//uint64_t start_part; // take from Index::index_parts_stats_vec
	bool load_for_search;
}; // ~class References
