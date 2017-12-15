#pragma once
/**
* FILE: read.hpp
* Created: Nov 06, 2017 Mon
* Wrapper of a Reads' record and its Match results
*/

#include <string>
#include <vector>
#include <algorithm> // std::find_if

#include "load_index.hpp"
#include "kvdb.hpp"
#include "traverse_bursttrie.hpp" // id_win
#include "ssw.hpp" // s_align2
#include "options.hpp"

class References; // forward

struct alignment_struct2
{
	uint32_t max_size; // max size of s_align array
	uint32_t size; // actual size of s_align array
	uint32_t min_index;
	uint32_t max_index;
	std::vector<s_align2> alignv;

	alignment_struct2() : max_size(0), size(0), min_index(0), max_index(0) {}

	alignment_struct2(std::string); // create from binary string

	alignment_struct2(uint32_t max_size, uint32_t size, uint32_t min, uint32_t max)
		: max_size(max_size), size(size), min_index(min), max_index(max) {}

	std::string toString(); // convert to binary string
	size_t getSize() { 
		size_t ret = sizeof(max_size) + sizeof(size) + sizeof(min_index) + sizeof(max_index);
		for (std::vector<s_align2>::iterator it = alignv.begin(); it != alignv.end(); ++it)
			ret += it->size();
		return ret;
	}
};

class Read {
public:
	int id = 0; // number of the read in the reads file
	bool isValid; // flags the record is not valid
	bool isEmpty; // flags the Read object is empty i.e. just a placeholder for copy assignment
	bool is03; // indicates Read::isequence is in 0..3 alphabet
	bool is04; // indicates Read:iseqeunce is in 0..4 alphabet. Seed search cannot proceed on 0-4 alphabet

	std::string header;
	std::string sequence;
	std::string quality; // "" (fasta) | "xxx..." (fastq)
	Format format; // fasta | fastq

	// calculated
	std::string isequence; // sequence in Integer alphabet: [A,C,G,T] -> [0,1,2,3]
	bool reversed = false;
	std::vector<int> ambiguous_nt; // positions of ambiguous nucleotides in the sequence (as defined in nt_table/load_index.cpp)

	// store in database ------------>
	int lastIndex; // last index number this read was aligned against. Set in Processor::callback
	int lastPart; // last part number this read was aligned against.  Set in Processor::callback
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
	std::vector<int8_t> scoring_matrix; // initScoringMatrix   orig: int8_t* scoring_matrix
	// <------------------------------ store in database

	const char complement[4] = { 3, 2, 1, 0 };

	Read()
		:
		isValid(false),
		isEmpty(true),
		is03(false),
		is04(false),
		lastIndex(-1),
		lastPart(-1)
	{
		if (num_alignments_gv > 0) num_alignments = num_alignments_gv;
		if (min_lis_gv > 0) best = min_lis_gv;
	}

	Read(int id, std::string header, std::string sequence, std::string quality, Format format)
		:
		id(id), header(std::move(header)), sequence(sequence),
		quality(quality), format(format), isEmpty(false)
	{
		validate();
		//seqToIntStr();
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
		is03 = that.is03;
		is04 = that.is04;
		header = that.header;
		sequence = that.sequence;
		quality = that.quality;
		format = that.format;
		isequence = that.isequence;
		reversed = that.reversed;
		ambiguous_nt = that.ambiguous_nt;
		lastIndex = that.lastIndex;
		lastPart = that.lastPart;
		hit = that.hit;
		hit_denovo = that.hit_denovo;
		null_align_output = that.null_align_output;
		max_SW_score = that.max_SW_score;
		num_alignments = that.num_alignments;
		readhit = that.readhit;
		best = that.best;
		id_win_hits = that.id_win_hits;
		hits_align_info = that.hits_align_info;
		scoring_matrix = that.scoring_matrix;
	}

	// copy assignment
	Read & operator=(const Read & that)
	{
		if (&that == this) return *this;

		//printf("Read copy assignment called\n");
		id = that.id;
		isValid = that.isValid;
		isEmpty = that.isEmpty;
		is03 = that.is03;
		is04 = that.is04;
		header = that.header;
		sequence = that.sequence;
		quality = that.quality;
		format = that.format;
		isequence = that.isequence;
		reversed = that.reversed;
		ambiguous_nt = that.ambiguous_nt;
		lastIndex = that.lastIndex;
		lastPart = that.lastPart;
		hit = that.hit;
		hit_denovo = that.hit_denovo;
		null_align_output = that.null_align_output;
		max_SW_score = that.max_SW_score;
		num_alignments = that.num_alignments;
		readhit = that.readhit;
		best = that.best;
		id_win_hits = that.id_win_hits;
		hits_align_info = that.hits_align_info;
		scoring_matrix = that.scoring_matrix;

		return *this; // by convention always return *this
	}

	void initScoringMatrix(Runopts & opts);

	// convert char "sequence" to 0..3 alphabet "isequence", and populate "ambiguous_nt"
	void seqToIntStr() 
	{
		for (std::string::iterator it = sequence.begin(); it != sequence.end(); ++it)
		{
			char c = (4 == nt_table[(int)*it]) ? 0 : nt_table[(int)*it];
			//isequence += nt_table[(int)*it];
			isequence += c;
			if (c == 0) { // ambiguous nt
				ambiguous_nt.push_back(static_cast<UINT>(isequence.size()) - 1); // i.e. add current position to the vector
			}
		}
		is03 = true;
	}

	// reverse complement the integer sequence
	void revIntStr() {
		std::reverse(isequence.begin(), isequence.end());
		for (int i = 0; i < isequence.length(); i++) {
			isequence[i] = complement[(int)isequence[i]]; // original: myread_rc[j] = complement[(int)*revcomp--]; paralleltraversal.cpp:975
			//isequence[i] = complement[isequence[i] - '0'];
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

	void init(Runopts & opts, KeyValueDatabase & kvdb)
	{
		validate();
		seqToIntStr();
		//unmarshallJson(kvdb); // get matches from Key-value database
		restoreFromDb(kvdb); // get matches from Key-value database
		initScoringMatrix(opts);
	}

	std::string matchesToJson(); // convert to Json string to store in DB

	std::string toString(); // convert to binary string to store in DB

	  // deserialize matches from string
	void restoreFromDb(KeyValueDatabase & kvdb);

	// deserialize matches from JSON and populate the read
	void unmarshallJson(KeyValueDatabase & kvdb);

	// flip isequence between 03 - 04 alphabets
	void flip34(Runopts & opts)
	{
		int val = is03 ? 4 : 0;
		if (ambiguous_nt.size() > 0) 
		{
			if (opts.forward)
			{
				for (uint32_t p = 0; p < ambiguous_nt.size(); p++)
				{
					isequence[ambiguous_nt[p]] = val;
				}
			}
			else
			{
				for (uint32_t p = 0; p < ambiguous_nt.size(); p++)
				{
					isequence[(isequence.length() - ambiguous_nt[p]) - 1] = val;
				}
			}
		}
	} // ~flip34

	void calcMismatchGapId(References & refs, int alignIdx, uint32_t & mismatches, uint32_t & gaps, uint32_t & id);
	std::string getSeqId() {
		std::string id = header.substr(0, header.find(' '));
		id.erase(id.begin(), std::find_if(id.begin(), id.end(), [](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));
		return id;
	}
}; // ~class Read
