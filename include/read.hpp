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
			   Laurent No�      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mika�l Salson    mikael.salson@lifl.fr
			   H�l�ne Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

/*
 * FILE: read.hpp
 * Created: Nov 06, 2017 Mon
 * Wrapper of a Reads' record and its Match results
 */

#pragma once

#include <string>
#include <sstream>
#include <vector>
#include <algorithm> // std::find_if
#include <tuple>

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
	std::vector<s_align2> alignv; // read alignments i.e. passing SW threshold

	// constructors
	alignment_struct2();
	alignment_struct2(std::string); // create from binary string

	// member functions
	std::string toString(); // convert to binary string
	size_t getSize();
	void clear();
};

/* 
 * 1. id_win_hits  std::vector<id_win> id_win_hits
 *	  array of positions of window hits on the reference sequence in given index/part.
 *	  Only used during alignment on a particular index/part. No need to store on disk.
 *	  Reset on each new index part
 *    [0] : {id = 568 win = 0 	}
 *	  ...
 *	  [4] : {id = 1248788 win = 72 }
 *	          |            |_k-mer start position on read
 *	          |_k-mer id on reference (index into 'positions_tbl')
 */
class Read 
{
public:
	std::string id; // Read ID: combination 'readsfile-number_read-number' e.g. 0_0, 1_103, etc.
	std::size_t read_num; // Read number in the reads file starting from 0
	std::size_t readfile_idx; // index into Runopts::readfiles
	bool isValid; // flags the record valid/non-valid
	bool isEmpty; // flags the Read object is empty i.e. just a placeholder for copy assignment
	bool is_too_short; // calculate for each reference file
	bool is03; // indicates Read::isequence is in 0..3 alphabet
	bool is04; // indicates Read:iseqeunce is in 0..4 alphabet. Seed search cannot proceed on 0-4 alphabet
	bool isRestored; // flags the read is restored from Database. See 'Read::restoreFromDb'
	BIO_FORMAT format; // fasta | fastq
	// store in database ------------>
	unsigned lastIndex; // last index number this read was aligned against. Set in Processor::callback
	unsigned lastPart; // last part number this read was aligned against.  Set in Processor::callback
	// matching results
	unsigned c_yid_ycov; // count of alignments passing both ID + COV
	unsigned n_yid_ncov; // count of alignments ID + !COV
	unsigned n_nid_ycov; // count of alignments !ID + COV
	unsigned n_denovo; // count of alignment !ID + !COV
	bool reversed; // indicates the read is reverse-complement i.e. 'revIntStr' was applied
	bool is_done; // all alignments have been found => stop searching
	bool is_hit; // true if at least one alignment 'SW_score >= min_SW_score' has been found
	bool is_new_hit; // indicates a new hit was found so the read has to be stored. Init to False before each index search. NO DB store.
	bool null_align_output; // flags NULL alignment was output to file (needs to be done once only)
	uint16_t max_SW_count; // count of alignments that have Max possible Smith-Waterman score for this read
	int32_t num_alignments; // counter of alignments to keep for reporting
	/*
	* count of read's k-mers that have DB matches. Used to filter the candidate references
	* that have a minimum required number of matching seeds (see opts.num_seeds)
	* 
	* TODO: It is incremented each time matches for a k-mer are found without regards to 
	* what references are involved. Whence it is possible for N k-mers to be matched each 
	* to a different reference, so that the 'hit_seeds' is N, but in reality each reference 
	* has only a single matching k-mer.
	* Not a problem because the matches on the refs are always calculated prior alignment.
	* This var is stil useful to ensure there are matches prior alignment.
	*/
	uint32_t hit_seeds;
	int32_t best; // init with opts.min_lis, see 'this.init'. NO DB store/restore (bug 51).
	std::string header;
	std::string sequence;
	std::string quality; // "" (fasta) | "xxx..." (fastq)
	// calculated
	std::string isequence; // sequence in Integer alphabet: [A,C,G,T] -> [0,1,2,3]
	std::vector<int> ambiguous_nt; // positions of ambiguous nucleotides in the sequence (as defined in nt_table/load_index.cpp)
	std::vector<id_win> id_win_hits; // [1] positions of kmer hits on the reference sequence in given index/part. NO DB store.
	alignment_struct2 alignment; // store in DB
	std::vector<int8_t> scoring_matrix; // initScoringMatrix   orig: int8_t* scoring_matrix  No DB store
	// <---- END store in database

public:
	Read();
	Read(std::string& readstr);
	Read(std::string id, std::size_t read_num);
	Read(std::string id, std::string header, std::string sequence, std::string quality, BIO_FORMAT format);
	Read(const Read & that); // copy constructor
	Read & operator=(const Read & that); // copy assignment
	//~Read();

public:
	void generate_id();
	void initScoringMatrix(int8_t match, int8_t mismatch, int8_t score_N);
	void validate();
	void clear();
	void init(Runopts& opts);
	/* convert to Json string to store in DB */
	// std::string matchesToJson();
	void unmarshallJson(KeyValueDatabase& kvdb);
	/* serialize to binary string to store in DB */
	std::string toBinString(); 
	bool load_db(KeyValueDatabase& kvdb);
	void seqToIntStr();
	void revIntStr();
	/* convert isequence to alphabetic form i.e. to A,C,G,T,N */
	std::string get04alphaSeq();
	/* flip isequence between 03 - 04 alphabets */
	void flip34();

	/*
	* count mismatches, gaps, matches, and calculate %ID, %COV given an alignment
	* @param  refs   references
	* @param  align  alignment to evaluate
	* @return tuple<mismatches, gaps, matches, %ID, %COV>
	*/
	std::tuple<uint32_t, uint32_t, uint32_t, double, double> calc_miss_gap_match(const References& refs, const s_align2& align);

	std::string getSeqId();
	uint32_t hashKmer(uint32_t pos, uint32_t len);
	bool from_string(std::string& readstr);
}; // ~class Read
