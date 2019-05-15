#pragma once
/**
* FILE: read.hpp
* Created: Nov 06, 2017 Mon
* Wrapper of a Reads' record and its Match results
*/

#include <string>
#include <sstream>
#include <vector>
#include <algorithm> // std::find_if

#include "kvdb.hpp"
#include "traverse_bursttrie.hpp" // id_win
#include "ssw.hpp" // s_align2
#include "options.hpp"

class References; // forward

struct alignment_struct2
{
	uint32_t max_size; // max size of alignv i.e. max number of alignments to store (see options '-N', '--best N') TODO: remove?
	uint32_t min_index; // index into alignv for the reference with lowest SW alignment score (see 'compute_lis_alignment')
	uint32_t max_index; // index into alignv for the reference with highest SW alignment score (see 'compute_lis_alignment')
	std::vector<s_align2> alignv;

	alignment_struct2() : max_size(0), min_index(0), max_index(0) {}

	alignment_struct2(std::string); // create from binary string

	std::string toString(); // convert to binary string
	size_t getSize() { 
		size_t ret = sizeof(min_index) + sizeof(max_index);
		for (std::vector<s_align2>::iterator it = alignv.begin(); it != alignv.end(); ++it)
			ret += it->size();
		return ret;
	}

	void clear()
	{
		max_size = 0;
		min_index = 0;
		max_index = 0;
		alignv.clear();
	}
};

class Read {
public:
	std::size_t id; // Read ID: hash combinations of read_num and readsfile
	std::size_t read_num; // Read number in the reads file. Use as key into Key-value Database.
	std::string readsfile; // reads file name
	bool isValid; // flags the record valid/non-valid
	bool isEmpty; // flags the Read object is empty i.e. just a placeholder for copy assignment
	bool is03; // indicates Read::isequence is in 0..3 alphabet
	bool is04; // indicates Read:iseqeunce is in 0..4 alphabet. Seed search cannot proceed on 0-4 alphabet
	bool isRestored; // flags the read is restored from Database. See 'Read::restoreFromDb'

	std::string header;
	std::string sequence;
	std::string quality; // "" (fasta) | "xxx..." (fastq)
	Format format; // fasta | fastq

	// calculated
	std::string isequence; // sequence in Integer alphabet: [A,C,G,T] -> [0,1,2,3]
	bool reversed = false; // indicates the read is reverse-complement i.e. 'revIntStr' was applied
	std::vector<int> ambiguous_nt; // positions of ambiguous nucleotides in the sequence (as defined in nt_table/load_index.cpp)

	// store in database ------------>
	unsigned int lastIndex; // last index number this read was aligned against. Set in Processor::callback
	unsigned int lastPart; // last part number this read was aligned against.  Set in Processor::callback
	// matching results
	bool hit = false; // indicates a match for this Read has been found
	bool hit_denovo = true; // hit & !(%Cov & %ID) TODO: change this to 'hit_cov_id' because it's set to true if !(%Cov & %ID) regardless of 'hit'
	bool null_align_output = false; // flags NULL alignment was output to file (needs to be done once only)
	uint16_t max_SW_count = 0; // count of matches that have Max Smith-Waterman score for this read
	int32_t num_alignments = 0; // number of alignments to output per read
	uint32_t readhit = 0; // number of seeds matches between read and database. (? readhit == id_win_hits.size)
	int32_t best = 0; // init with opts.min_lis, see 'this.init'. Don't DB store/restore (bug 51).

	// array of positions of window hits on the reference sequence in given index/part. 
	// Only used during alignment on a particular index/part. No need to store on disk.
	// Reset on each new index part
	// [0] : {id = 568 win = 0 	}
	// ...
	// [4] : {id = 1248788 win = 72 }
	//        |            |_k-mer start position on read
	//        |_k-mer id on reference (index into 'positions_tbl')
	std::vector<id_win> id_win_hits;

	alignment_struct2 hits_align_info; // stored in DB

	std::vector<int8_t> scoring_matrix; // initScoringMatrix   orig: int8_t* scoring_matrix
	// <------------------------------ store in database

	Read();
	Read(int id, std::string header, std::string sequence, std::string quality, Format format);
	~Read();
	Read(const Read & that); // copy constructor
	Read & operator=(const Read & that); // copy assignment

	void generate_id();
	void initScoringMatrix(long match, long mismatch, long score_N);

	void validate();

	void clear();

	void init(Runopts & opts);

	std::string matchesToJson(); // convert to Json string to store in DB

	void unmarshallJson(KeyValueDatabase & kvdb); // deserialize matches from JSON and populate the read

	std::string toString(); // convert to binary string to store in DB

	bool load_db(KeyValueDatabase & kvdb); // deserialize matches from string

	void seqToIntStr();

	void revIntStr(); // reverse complement the integer sequence in 03 encoding

	std::string get04alphaSeq(); // convert isequence to alphabetic form i.e. to A,C,G,T,N

	void flip34(); // flip isequence between 03 - 04 alphabets

	void calcMismatchGapId(References & refs, int alignIdx, uint32_t & mismatches, uint32_t & gaps, uint32_t & id);

	std::string getSeqId();

	uint32_t hashKmer(uint32_t pos, uint32_t len);
}; // ~class Read
