/**
 * FILE: read.cpp
 * Created: Nov 26, 2017 Sun
 * @copyright 2016-20 Clarity Genomics BVBA
 */
#include <filesystem>

// 3rd party
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

// SMR
#include "read.hpp"
#include "references.hpp"

alignment_struct2::alignment_struct2() : max_size(0), min_index(0), max_index(0) 
{}

/**
 * construct from the binary string stored in DB
 */
alignment_struct2::alignment_struct2(std::string bstr)
{
	size_t offset = 0;
	// min_index
	std::memcpy(static_cast<void*>(&min_index), bstr.data() + offset, sizeof(min_index));
	offset += sizeof(min_index);
	// max_index
	std::memcpy(static_cast<void*>(&max_index), bstr.data() + offset, sizeof(max_index));
	offset += sizeof(max_index);

	// alignv vector<s_align2>
	size_t alignv_size = 0; // size of alignv vector
	std::memcpy(static_cast<void*>(&alignv_size), bstr.data() + offset, sizeof(alignv_size));
	offset += sizeof(alignv_size);

	for (int i = 0; i < alignv_size; i++)
	{
		size_t s_align_size = 0; // size of a_align struct
		std::memcpy(static_cast<void*>(&s_align_size), bstr.data() + offset, sizeof(s_align_size));
		offset += sizeof(s_align_size);
		std::string s_align_str(bstr.data() + offset, bstr.data() + offset + s_align_size); // part of the string encoding the current s_align struct
		offset += s_align_size;
		alignv.push_back(s_align2(s_align_str));
	}
}

std::string alignment_struct2::toString()
{
	std::string buf;
	size_t offset = 0;

	std::copy_n(static_cast<char*>(static_cast<void*>(&min_index)), sizeof(min_index), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&max_index)), sizeof(max_index), std::back_inserter(buf));

	// alignv
	size_t vsize = alignv.size();
	std::copy_n(static_cast<char*>(static_cast<void*>(&vsize)), sizeof(vsize), std::back_inserter(buf)); // add vector size

	// alignv
	std::string s_align2_str;
	size_t s_align2_strlen = 0;
	for (auto it = alignv.begin(); it < alignv.end(); ++it) {
		s_align2_str = it->toString(); // get string
		s_align2_strlen = s_align2_str.size(); // get string size
		std::copy_n(static_cast<char*>(static_cast<void*>(&s_align2_strlen)), sizeof(s_align2_strlen), std::back_inserter(buf)); // add size to buffer
		buf += s_align2_str; // add string to buffer
	}

	return buf;
} // ~toString

size_t alignment_struct2::getSize() {
	size_t ret = sizeof(min_index) + sizeof(max_index);
	for (std::vector<s_align2>::iterator it = alignv.begin(); it != alignv.end(); ++it)
		ret += it->size();
	return ret;
}

void alignment_struct2::clear()
{
	max_size = 0;
	min_index = 0;
	max_index = 0;
	alignv.clear();
}

Read::Read()
	:
	id(""),
	read_num(0),
	readfile_idx(0),
	isValid(false),
	isEmpty(true),
	is03(false),
	is04(false),
	isRestored(false),
	lastIndex(0),
	lastPart(0),
	reversed(false),
	is_hit(false),
	is_denovo(true),
	null_align_output(false),
	max_SW_count(0),
	num_alignments(0),
	readhit(0),
	best(0),
	format(Format::FASTA)
{}

Read::Read(std::string& readstr)
	:
	Read()
{
	isEmpty = !from_string(readstr);
}

Read::Read(std::string id, std::string& read)
	:
	Read()
	{ id = id; }

Read::Read(std::string id, std::string header, std::string sequence, std::string quality, Format format)
	:
	id(id), header(std::move(header)), sequence(sequence),
	quality(quality), format(format), isEmpty(false)
{
	validate();
}

Read::~Read() {}

// copy constructor
Read::Read(const Read & that)
{
	id = that.id;
	read_num = that.read_num;
	readfile_idx = that.readfile_idx;
	isValid = that.isValid;
	isEmpty = that.isEmpty;
	is03 = that.is03;
	is04 = that.is04;
	isRestored = that.isRestored;
	header = that.header;
	sequence = that.sequence;
	quality = that.quality;
	format = that.format;
	isequence = that.isequence;
	reversed = that.reversed;
	ambiguous_nt = that.ambiguous_nt;
	lastIndex = that.lastIndex;
	lastPart = that.lastPart;
	is_hit = that.is_hit;
	is_denovo = that.is_denovo;
	null_align_output = that.null_align_output;
	max_SW_count = that.max_SW_count;
	num_alignments = that.num_alignments;
	readhit = that.readhit;
	best = that.best;
	id_win_hits = that.id_win_hits;
	alignment = that.alignment;
	scoring_matrix = that.scoring_matrix;
}

// copy assignment
Read & Read::operator=(const Read & that)
{
	if (&that == this) return *this; // else *this = that

	//printf("Read copy assignment called\n");
	id = that.id;
	read_num = that.read_num;
	readfile_idx = that.readfile_idx;
	isValid = that.isValid;
	isEmpty = that.isEmpty;
	is03 = that.is03;
	is04 = that.is04;
	isRestored = that.isRestored;
	header = that.header;
	sequence = that.sequence;
	quality = that.quality;
	format = that.format;
	isequence = that.isequence;
	reversed = that.reversed;
	ambiguous_nt = that.ambiguous_nt;
	lastIndex = that.lastIndex;
	lastPart = that.lastPart;
	is_hit = that.is_hit;
	is_denovo = that.is_denovo;
	null_align_output = that.null_align_output;
	max_SW_count = that.max_SW_count;
	num_alignments = that.num_alignments;
	readhit = that.readhit;
	best = that.best;
	id_win_hits = that.id_win_hits;
	alignment = that.alignment;
	scoring_matrix = that.scoring_matrix;

	return *this; // by convention always return *this
} // ~Read::operator=

/** 
 * Generate ID of the read
 */
void Read::generate_id()
{
	std::stringstream ss;
	id.clear();
	ss << readfile_idx << "_" << read_num; // << std::filesystem::path(opts.readfiles[readfile_num]).filename();
	char buf[4096];
	while (ss.read(buf, sizeof(buf)))
		id.append(buf, sizeof(buf));
	id.append(buf, ss.gcount());
	//std::hash<std::string> hash_fn;
	//id = hash_fn(ss.str());
} // ~Read::generate_id

/**
 * 5 options are used here, which would make this method to take 7 args => use Runopts as arg
 */
void Read::init(Runopts & opts)
{
	if (opts.num_alignments > 0) this->num_alignments = opts.num_alignments;
	if (opts.min_lis > 0) this->best = opts.min_lis;
	validate();
	seqToIntStr();
	initScoringMatrix(opts.match, opts.mismatch, opts.score_N);
} // ~Read::init

// initialize Smith-Waterman scoring matrix for genome sequences
void Read::initScoringMatrix(int8_t match, int8_t mismatch, int8_t score_N)
{
	int8_t val = 0;
	for (int l = 0; l < 4; ++l)
	{
		for (int m = 0; m < 4; ++m) {
			val = l == m ? match : mismatch; // weight_match : weight_mismatch (must be negative)
			scoring_matrix.push_back(val);
		}
		scoring_matrix.push_back(score_N); // ambiguous base
	}
	for (int m = 0; m < 5; ++m) {
		scoring_matrix.push_back(score_N); // ambiguous base
	}
}

void Read::validate() {
	std::stringstream ss;
	if (sequence.size() > MAX_READ_LEN)
	{
		ERR("Read ID: ", id, " Header: ", header, " Sequence length: ", sequence.size(), " > ", 
			MAX_READ_LEN, " nt \n", "  Please check your reads or contact the authors.");
		exit(EXIT_FAILURE);
	}
	isValid = true;
} // ~Read::validate

void Read::clear()
{
	id.clear();
	isValid = false;
	isEmpty = true;
	is03 = false;
	is04 = false;
	header.clear();
	sequence.clear();
	quality.clear();
	isequence.clear();
	reversed = false;
	ambiguous_nt.clear();
	isRestored = false;
	lastIndex = 0;
	lastPart = 0;
	is_hit = false;
	is_denovo = true;
	null_align_output = false;
	max_SW_count = 0;
	num_alignments = 0;
	readhit = 0;
	best = 0;
	id_win_hits.clear();
	alignment.clear();
	scoring_matrix.clear();
} // ~Read::clear

// convert char "sequence" to 0..3 alphabet "isequence", and populate "ambiguous_nt"
void Read::seqToIntStr()
{
	for (std::string::iterator it = sequence.begin(); it != sequence.end(); ++it)
	{
		char c = nt_table[(int)*it];
		if (c == 4) // ambiguous nt. 4 is max value in nt_table
		{
			ambiguous_nt.push_back(static_cast<int>(isequence.size())); // i.e. add current position to the vector
			c = 0;
		}
		isequence += c;
	}
	is03 = true;
}

/* reverse complement the integer sequence in 03 encoding */
void Read::revIntStr() 
{
	std::reverse(isequence.begin(), isequence.end());
	for (int i = 0; i < isequence.length(); i++) {
		isequence[i] = complement[(int)isequence[i]];
	}
	reversed = !reversed;
}

// convert isequence to alphabetic form i.e. to A,C,G,T,N
std::string Read::get04alphaSeq() {
	//bool rev03 = false; // mark whether to revert back to 03
	std::string seq;
	if (is03) flip34();
	// convert to alphabetic
	for (int i = 0; i < isequence.size(); ++i)
		seq += nt_map[(int)isequence[i]];

	//if (rev03) flip34();
	return seq;
}

std::string Read::getSeqId() {
	// part of the header from start till first space.
	std::string id = header.substr(0, header.find(' '));
	// remove '>' or '@'
	id.erase(id.begin(), std::find_if(id.begin(), id.end(), [](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));
	return id;
}

// flip isequence between 03 - 04 alphabets
void Read::flip34()
{
	if (ambiguous_nt.size() > 0)
	{
		int val = is03 ? 4 : 0;
		if (reversed)
		{
			for (uint32_t p = 0; p < ambiguous_nt.size(); p++)
			{
				isequence[(isequence.length() - ambiguous_nt[p]) - 1] = val;
			}
		}
		else
		{
			for (uint32_t p = 0; p < ambiguous_nt.size(); p++)
			{
				isequence[ambiguous_nt[p]] = val;
			}
		}
		is03 = !is03;
		is04 = !is04;
	}
} // ~flip34

/* convert to Json string to store in DB */
std::string Read::matchesToJson() {
	rapidjson::StringBuffer sbuf;
	rapidjson::Writer<rapidjson::StringBuffer> writer(sbuf);

	writer.StartObject();
	writer.Key("hit");
	writer.Bool(is_hit);
	writer.Key("is_denovo");
	writer.Bool(is_denovo);
	writer.Key("null_align_output");
	writer.Bool(null_align_output);
	writer.Key("max_SW_count");
	writer.Uint(max_SW_count);
	writer.Key("num_alignments");
	writer.Int(num_alignments);

	writer.Key("alignment");
	writer.StartObject();
	writer.String("max_size");
	writer.Uint(10);
	writer.EndObject();

	writer.EndObject();

	return sbuf.GetString();
} // ~Read::matchesToJsonString

/* 
 * serialize to binary string to store in DB 
 */
std::string Read::toBinString()
{
	if (alignment.alignv.size() == 0)
		return "";

	// hit, hit_denovo, null_align_output, max_SW_count, num_alignments, readhit, best
	std::string buf;
	std::copy_n(static_cast<char*>(static_cast<void*>(&lastIndex)), sizeof(lastIndex), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&lastPart)), sizeof(lastPart), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&is_hit)), sizeof(is_hit), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&is_denovo)), sizeof(is_denovo), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&null_align_output)), sizeof(null_align_output), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&max_SW_count)), sizeof(max_SW_count), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&num_alignments)), sizeof(num_alignments), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&readhit)), sizeof(readhit), std::back_inserter(buf));
	//std::copy_n(static_cast<char*>(static_cast<void*>(&best)), sizeof(best), std::back_inserter(buf));

	// id_win_hits vector - TODO: remove?
#if 0
	size_t id_win_hits_size = id_win_hits.size();
	std::copy_n(static_cast<char*>(static_cast<void*>(&id_win_hits_size)), sizeof(id_win_hits_size), std::back_inserter(buf)); // add vector size
	for (auto it = id_win_hits.begin(); it != id_win_hits.end(); ++it) buf += it->toString(); // add values
#endif
	// alignment
	std::string alignment_str = alignment.toString();
	size_t alignment_size = alignment_str.length();
	std::copy_n(static_cast<char*>(static_cast<void*>(&alignment_size)), sizeof(alignment_size), std::back_inserter(buf)); // add size
	buf += alignment_str; //  add string

	return buf;
} // ~Read::toBinString

/* 
 * load read alignment data from DB 
 */
bool Read::load_db(KeyValueDatabase & kvdb)
{
	int id_win_hits_len = 0;
	std::string bstr = kvdb.get(id);
	if (bstr.size() == 0) { isRestored = false; return isRestored; }
	size_t offset = 0;

	std::memcpy(static_cast<void*>(&lastIndex), bstr.data() + offset, sizeof(lastIndex));
	offset += sizeof(lastIndex);

	std::memcpy(static_cast<void*>(&lastPart), bstr.data() + offset, sizeof(lastPart));
	offset += sizeof(lastPart);

	std::memcpy(static_cast<void*>(&is_hit), bstr.data() + offset, sizeof(is_hit));
	offset += sizeof(is_hit);

	std::memcpy(static_cast<void*>(&is_denovo), bstr.data() + offset, sizeof(is_denovo));
	offset += sizeof(is_denovo);

	std::memcpy(static_cast<void*>(&null_align_output), bstr.data() + offset, sizeof(null_align_output));
	offset += sizeof(null_align_output);

	std::memcpy(static_cast<void*>(&max_SW_count), bstr.data() + offset, sizeof(max_SW_count));
	offset += sizeof(max_SW_count);

	std::memcpy(static_cast<void*>(&num_alignments), bstr.data() + offset, sizeof(num_alignments));
	offset += sizeof(num_alignments);

	std::memcpy(static_cast<void*>(&readhit), bstr.data() + offset, sizeof(readhit));
	offset += sizeof(readhit);

	//std::memcpy(static_cast<void*>(&best), bstr.data() + offset, sizeof(best));
	//offset += sizeof(best);

	// std::vector<id_win> id_win_hits TODO: remove?
#if 0
	size_t id_win_hits_size = 0;
	std::memcpy(static_cast<void*>(&id_win_hits_size), bstr.data() + offset, sizeof(id_win_hits_size));
	offset += sizeof(id_win_hits_size);
	std::string id_win_str;
	for (int i = 0; i < id_win_hits_size; ++i)
	{
		//id_win_str.assign(sizeof(id_win::id) + sizeof(id_win::win), 0);
		//std::memcpy(static_cast<void*>(&id_win_str[0]), bstr.data() + offset, sizeof(id_win::id) + sizeof(id_win::win));
		std::copy_n(bstr.data() + offset, sizeof(id_win::id) + sizeof(id_win::win), std::back_inserter(id_win_str));
		id_win_hits.push_back(id_win(id_win_str));
		offset += (sizeof(id_win::id) + sizeof(id_win::win));
		id_win_str.clear();
	}
#endif
	// alignment_struct2 alignment
	size_t alignment_size = 0;
	std::memcpy(static_cast<void*>(&alignment_size), bstr.data() + offset, sizeof(alignment_size));
	offset += sizeof(alignment_size);
	std::string alignment_str(bstr.data() + offset, bstr.data() + offset + alignment_size);
	alignment_struct2 alignstruct(alignment_str);
	alignment = alignstruct;
	offset += alignment_size;

	isRestored = true;
	return isRestored;
} // ~Read::load_db

/* deserialize matches from JSON and populate the read */
void Read::unmarshallJson(KeyValueDatabase & kvdb)
{
	INFO("Read::unmarshallJson: Not yet Implemented");
}



/**
* Prototype: paralleltraversal lines 1531..1555
* Calculate Mismatches, Gaps, and ID
*
* @param IN Refs  references
* @param IN Read
* @param IN alignIdx index into Read.alignment.alignv
* @param OUT mismatches  calculated here for the given Read Alignment
* @param OUT gaps
* @param OUT id   matched characters
*/
void Read::calcMismatchGapId(References & refs, int alignIdx, uint32_t & mismatches, uint32_t & gaps, uint32_t & id)
{
	if (alignIdx >= alignment.alignv.size()) return; // index exceeds the size of the alignment vector

	mismatches = 0; // count of mismatched characters
	gaps = 0; // count of gaps
	id = 0; // count of matched characters

	int32_t qb = alignment.alignv[alignIdx].ref_begin1; // index of the first char in the reference matched part
	int32_t pb = alignment.alignv[alignIdx].read_begin1; // index of the first char in the read matched part

	std::string refseq = refs.buffer[alignment.alignv[alignIdx].ref_seq].sequence;

	for (uint32_t cidx = 0; cidx < alignment.alignv[alignIdx].cigar.size(); ++cidx)
	{
		uint32_t letter = 0xf & alignment.alignv[alignIdx].cigar[cidx]; // 4 low bits
		uint32_t length = (0xfffffff0 & alignment.alignv[alignIdx].cigar[cidx]) >> 4; // high 28 bits i.e. 32-4=28
		if (letter == 0)
		{
			for (uint32_t u = 0; u < length; ++u)
			{
				if (refseq[qb] != isequence[pb]) ++mismatches;
				else ++id;
				++qb;
				++pb;
			}
		}
		else if (letter == 1)
		{
			pb += length;
			gaps += length;
		}
		else
		{
			qb += length;
			gaps += length;
		}
	}
} // ~Read::calcMismatchGapId

/* 
 * Calculate the numerical value (hash) of a kmer given its position on the read and its length.
 * The hash is just a numeric value formed by the chars of a string consisting of '0','1','2','3'
 * e.g. "2233012" -> b10.1011.1100.0110 = x2BC6 = 11206
 *
 * Use to lookup the 'Index::lookup_tbl'
 *
 * @param pos  Kmer position on the read
 * @param len  Kmer Length
 */
uint32_t Read::hashKmer(uint32_t pos, uint32_t len)
{
	uint32_t hash = 0;
	char *pKmer = &isequence[pos];
	for (uint32_t i = 0; i < len; i++)
	{
		(hash <<= 2) |= (uint32_t)(*pKmer);
		++pKmer;
	}
	return hash;
}

/* 
 * @param readstr 'read_id \n header \n sequence [\n quality]'
 */
bool Read::from_string(std::string& readstr)
{
	bool is_ok = true;;
	std::stringstream ss(readstr);
	std::string line;
	for (int i = 0; std::getline(ss, line, '\n'); ++i) {
		if (i == 0) {
			id = line;
			auto pos = line.find_first_of('_');
			readfile_idx = std::atoi(line.substr(0, pos+1).data());
			read_num = std::atoi(line.substr(pos+1, line.size() - pos -1).data());
		}
		else if (i == 1) {
			format = line.front() == FASTA_HEADER_START ? Format::FASTA : Format::FASTQ;
			header = line;
		}
		else if (i == 2) {
			sequence = line;
		}
		else if (i == 3) {
			if (format == Format::FASTQ) {
				quality = line;
			}
			else {
				ERR("unexpected number of lines in fasta read: ", readstr);
				is_ok = false;
			}
		}
		else {
			ERR("unexpected number of lines in read: ", readstr);
			is_ok = false;
		}
	}
	return is_ok;
} // ~Read::from_string