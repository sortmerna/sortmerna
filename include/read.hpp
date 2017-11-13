#pragma once
/**
* FILE: read.hpp
* Created: Nov 06, 2017 Mon
* Wrapper of a Reads' record and its Match results
*/

#include <string>
#include <vector>

#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

#include "load_index.hpp"
#include "kvdb.hpp"


struct alignment_struct2
{
	uint32_t max_size; // max size of s_align array
	uint32_t size; // actual size of s_align array
	uint32_t min_index;
	uint32_t max_index;
	std::vector<s_align2> alignv;

	alignment_struct2() : max_size(0), size(0), min_index(0), max_index(0) {}
	alignment_struct2(uint32_t max_size, uint32_t size, uint32_t min, uint32_t max)
		: max_size(max_size), size(size), min_index(min), max_index(max) {}

	// copy constructor
	alignment_struct2(const alignment_struct2 & that)
	{}

	// copy assignment
	alignment_struct2 & operator=(const alignment_struct2 & that)
	{
		if (this != &that)
		{
		}
		return *this;
	}

	// convert to binary string
	std::string toString()
	{
		int bufsize = 4 * sizeof(max_size);
		std::string buf(bufsize, 0);
		int bufidx = 0;
		// max_size, size, min_index, max_index
		char * pch = reinterpret_cast<char *>(&max_size);
		for (int i = 0; i < 4 * sizeof(max_size); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// alignv
		for (auto it = alignv.begin(); it < alignv.end(); ++it) buf.append(it->toString());

		return buf;
	}
};

class Read {
public:
	int id = 0; // number of the read in the reads file
	bool isValid; // flags the record is not valid
	bool isEmpty; // flags the Read object is empty i.e. just a placeholder for copy assignment

	std::string header;
	std::string sequence;
	std::string quality; // "" (fasta) | "xxx..." (fastq)
	Format format; // fasta | fastq

	// calculated
	std::string seq_int_str; // sequence in Integer alphabet: [A,C,G,T] -> [0,1,2,3]
	bool reversed = false;
	std::vector<int> ambiguous_nt; // positions of ambiguous nucleotides in the sequence (as defined in nt_table/load_index.cpp)

	// store in database ------------>
	// matching results
	bool hit = false; // indicates that a match for this Read has been found
	bool hit_denovo = true;
	bool null_align_output = false; // flags NULL alignment was output to file (needs to be done once only)
	uint16_t max_SW_score = 0; // Max Smith-Waterman score
	int32_t num_alignments = 0; // number of alignments to output per read
	uint32_t readhit = 0; // number of seeds matches between read and database. Total number of hits?
	int32_t best = 0; // init with min_lis_gv

	// need custom destructor, copy constructor, and copy assignment
	std::vector<id_win> id_win_hits; // array of positions of window hits on the reference sequence
	alignment_struct2 hits_align_info;
	std::vector<int8_t> scoring_matrix;
	//int8_t* scoring_matrix = (int8_t*)calloc(25, sizeof(int8_t));
	//int8_t* ss = new int8_t[25];
	//std::unique_ptr<int8_t[]> scoring_matrix2(new int8_t[25]);
	// <------------------------------ store in database

	const char complement[4] = { '3','2','1','0' };

	Read() : isValid(false), isEmpty(true), scoring_matrix(25, 0)
	{
		if (num_alignments_gv > 0) num_alignments = num_alignments_gv;
		if (min_lis_gv > 0) best = min_lis_gv;
		// create new instance of alignments
		//hits_align_info.max_size = 0;
		//hits_align_info.size = 0;
		//hits_align_info.min_index = 0;
		//hits_align_info.max_index = 0;
		//hits_align_info.ptr = new s_align[1](); // see alignment.cpp
		//hits_align_info.ptr->cigar = 0;
		//hits_align_info.ptr->cigar = new uint32_t[1];
		//hits_align_info.ptr->cigarLen = 0;
		//hits_align_info.ptr->index_num = 0;
		//hits_align_info.ptr->part = 0;
		//hits_align_info.ptr->readlen = 0;
		//hits_align_info.ptr->read_begin1 = 0;
		//hits_align_info.ptr->read_end1 = 0;
		//hits_align_info.ptr->ref_begin1 = 0;
		//hits_align_info.ptr->ref_end1 = 0;
		//hits_align_info.ptr->ref_seq = 0;
		//hits_align_info.ptr->score1 = 0;
		//hits_align_info.ptr->strand = 0;
	}

	Read(int id, std::string header, std::string sequence, std::string quality, Format format)
		:
		id(id), header(std::move(header)), sequence(sequence),
		quality(quality), format(format), isEmpty(false)
	{
		validate();
		seqToIntStr();
		//		initScoringMatrix(opts);
	}

	~Read() {
		//		if (hits_align_info.ptr != 0) {
		//	delete[] read.hits_align_info.ptr->cigar;
		//	delete read.hits_align_info.ptr;
		//			delete hits_align_info.ptr; 
		//		}
		//		free(scoring_matrix);
		//		scoring_matrix = 0;
	}

	// copy constructor
	Read(const Read & that)
	{
		id = that.id;
		isValid = that.isValid;
		isEmpty = that.isEmpty;
		header = that.header;
		sequence = that.sequence;
		seq_int_str = that.seq_int_str;
	}

	// copy assignment
	Read & operator=(const Read & that)
	{
		if (&that == this) return *this;

		printf("Read copy assignment called\n");
		id = that.id;
		isValid = that.isValid;
		isEmpty = that.isEmpty;
		header = that.header;
		sequence = that.sequence;
		seq_int_str = that.seq_int_str;

		return *this; // by convention always return *this
	}

	//	void initScoringMatrix(Runopts & opts);

	// convert sequence to "sequenceInt" and populate "ambiguous_nt"
	void seqToIntStr() {
		for (std::string::iterator it = sequence.begin(); it != sequence.end(); ++it)
		{
			char c = ('4' == nt_table[(int)*it]) ? '0' : nt_table[(int)*it];
			seq_int_str.append(1, nt_table[(int)*it]);
			if (c == '0') { // ambiguous nt
				ambiguous_nt.push_back(static_cast<UINT>(seq_int_str.size()) - 1); // i.e. add current position to the vector
			}
		}
	}

	// reverse complement the integer sequence
	void revIntStr() {
		std::reverse(seq_int_str.begin(), seq_int_str.end());
		for (int i = 0; i < seq_int_str.length(); i++) {
			seq_int_str[i] = complement[seq_int_str[i] - '0'];
		}
		reversed = true;
	}

	void validate() {
		if (sequence.size() > READLEN)
		{
			fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] at least one of your reads is > %d nt \n",
				startColor, "\033[0m", __LINE__, __FILE__, READLEN);
			fprintf(stderr, "  Please check your reads or contact the authors.\n");
			exit(EXIT_FAILURE);
		}
		isValid = true;
	} // ~validate

	void clear()
	{
		header.clear();
		sequence.clear();
		quality.clear();
		isValid = false;
		isEmpty = true;
	}

	void init(KeyValueDatabase & kvdb)
	{
		validate();
		seqToIntStr();
		unmarshallJson(kvdb); // get matches from Key-value database
	}


	std::string matchesToJson() {
		rapidjson::StringBuffer sbuf;
		rapidjson::Writer<rapidjson::StringBuffer> writer(sbuf);

		writer.StartObject();
		writer.Key("hit");
		writer.Bool(hit);
		writer.Key("hit_denovo");
		writer.Bool(hit_denovo);
		writer.Key("null_align_output");
		writer.Bool(null_align_output);
		writer.Key("max_SW_score");
		writer.Uint(max_SW_score);
		writer.Key("num_alignments");
		writer.Int(num_alignments);

		writer.Key("hits_align_info");
		writer.StartObject();
		writer.String("max_size");
		writer.Uint(10);
		writer.EndObject();

		writer.EndObject();

		return sbuf.GetString();
	} // ~Read::matchesToJsonString

	  // convert to binary string whatever needs to be stored in DB
	std::string toString() {
		if (hits_align_info.alignv.size() == 0)
			return "";

		// hit, hit_denovo, null_align_output, max_SW_score, num_alignments, readhit, best
		int bufsize = 3 * sizeof(hit) + sizeof(max_SW_score) + sizeof(num_alignments) + sizeof(readhit) + sizeof(best);
		std::string buf(bufsize, 0);
		int bufidx = 0;
		char* pch = reinterpret_cast<char *>(&hit);
		for (int i = 0; i < 4 * sizeof(bufsize); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// id_win_hits
		for (auto it = id_win_hits.begin(); it != id_win_hits.end(); ++it) buf.append(it->toString());
		// hits_align_info
		buf.append(hits_align_info.toString());
		// std::vector<int8_t> scoring_matrix;
		for (auto it = scoring_matrix.begin(); it != scoring_matrix.end(); ++it) buf.append(1, *it);

		return buf;
	} // ~Read::toString

	  // deserialize matches from string
	void unmarshallString(std::string matchStr);

	// deserialize matches from JSON and populate the read
	void unmarshallJson(KeyValueDatabase & kvdb);
}; // ~class Read
