/**
 * FILE: read.cpp
 * Created: Nov 26, 2017 Sun
 */

// 3rd party
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

// SMR
#include "read.hpp"
#include "references.hpp"

/**
 * construct from the binary string stored in DB
 */
alignment_struct2::alignment_struct2(std::string bstr)
{
	size_t offset = 0;
	// max_size
	std::memcpy(static_cast<void*>(&max_size), bstr.data(), sizeof(max_size));
	offset += sizeof(size);
	// size
	std::memcpy(static_cast<void*>(&size), bstr.data() + offset, sizeof(size));
	offset += sizeof(size);
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

	std::copy_n(static_cast<char*>(static_cast<void*>(&max_size)), sizeof(max_size), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&size)), sizeof(size), std::back_inserter(buf));
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

 // initialize Smith-Waterman scoring matrix for genome sequences
void Read::initScoringMatrix(Runopts & opts)
{
	int l, k, m;
	int8_t val = 0;
	for (l = k = 0; l < 4; ++l)
	{
		for (m = 0; m < 4; ++m) {
			val = l == m ? opts.match : opts.mismatch; // weight_match : weight_mismatch (must be negative)
			scoring_matrix.push_back(val);
		}
		scoring_matrix.push_back(opts.score_N); // ambiguous base
	}
	for (m = 0; m < 5; ++m) {
		scoring_matrix.push_back(opts.score_N); // ambiguous base
	}
}

std::string Read::matchesToJson() {
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

std::string Read::toString()
{
	if (hits_align_info.alignv.size() == 0)
		return "";

	// hit, hit_denovo, null_align_output, max_SW_score, num_alignments, readhit, best
	std::string buf;
	std::copy_n(static_cast<char*>(static_cast<void*>(&lastIndex)), sizeof(lastIndex), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&lastPart)), sizeof(lastPart), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&hit)), sizeof(hit), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&hit_denovo)), sizeof(hit_denovo), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&null_align_output)), sizeof(null_align_output), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&max_SW_score)), sizeof(max_SW_score), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&num_alignments)), sizeof(num_alignments), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&readhit)), sizeof(readhit), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&best)), sizeof(best), std::back_inserter(buf));

	// id_win_hits vector
	size_t id_win_hits_size = id_win_hits.size();
	std::copy_n(static_cast<char*>(static_cast<void*>(&id_win_hits_size)), sizeof(id_win_hits_size), std::back_inserter(buf)); // add vector size
	for (auto it = id_win_hits.begin(); it != id_win_hits.end(); ++it) buf += it->toString(); // add values

	// hits_align_info
	std::string hits_align_info_str = hits_align_info.toString();
	size_t hits_align_info_size = hits_align_info_str.length();
	std::copy_n(static_cast<char*>(static_cast<void*>(&hits_align_info_size)), sizeof(hits_align_info_size), std::back_inserter(buf)); // add size
	buf += hits_align_info_str; //  add string

	return buf;
} // ~Read::toString

void Read::restoreFromDb(KeyValueDatabase & kvdb)
{
	int id_win_hits_len = 0;
	std::string bstr = kvdb.get(std::to_string(id));
	if (bstr.size() == 0) return;
	size_t offset = 0;

	std::memcpy(static_cast<void*>(&lastIndex), bstr.data() + offset, sizeof(lastIndex));
	offset += sizeof(lastIndex);

	std::memcpy(static_cast<void*>(&lastPart), bstr.data() + offset, sizeof(lastPart));
	offset += sizeof(lastPart);

	std::memcpy(static_cast<void*>(&hit), bstr.data() + offset, sizeof(hit));
	offset += sizeof(hit);

	std::memcpy(static_cast<void*>(&hit_denovo), bstr.data() + offset, sizeof(hit_denovo));
	offset += sizeof(hit_denovo);

	std::memcpy(static_cast<void*>(&null_align_output), bstr.data() + offset, sizeof(null_align_output));
	offset += sizeof(null_align_output);

	std::memcpy(static_cast<void*>(&max_SW_score), bstr.data() + offset, sizeof(max_SW_score));
	offset += sizeof(max_SW_score);

	std::memcpy(static_cast<void*>(&num_alignments), bstr.data() + offset, sizeof(num_alignments));
	offset += sizeof(num_alignments);

	std::memcpy(static_cast<void*>(&readhit), bstr.data() + offset, sizeof(readhit));
	offset += sizeof(readhit);

	std::memcpy(static_cast<void*>(&best), bstr.data() + offset, sizeof(best));
	offset += sizeof(best);

	// std::vector<id_win> id_win_hits
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

	// alignment_struct2 hits_align_info
	size_t hits_align_info_size = 0;
	std::memcpy(static_cast<void*>(&hits_align_info_size), bstr.data() + offset, sizeof(hits_align_info_size));
	offset += sizeof(hits_align_info_size);
	std::string hits_align_info_str(bstr.data() + offset, bstr.data() + offset + hits_align_info_size);
	alignment_struct2 alignstruct(hits_align_info_str);
	hits_align_info = alignstruct;
	offset += hits_align_info_size;
} // ~Read::restoreFromDb

void Read::unmarshallJson(KeyValueDatabase & kvdb)
{
	printf("Read::unmarshallJson: Not yet Implemented\n");
}



/**
* Prototype: paralleltraversal lines 1531..1555
* Calculate Mismatches, Gaps, and ID
*
* @param IN Refs  references
* @param IN Read
* @param IN alignIdx index into Read.hits_align_info.alignv
* @param OUT mismatches  calculated here for the given Read Alignment
* @param OUT gaps
* @param OUT id
*/
void Read::calcMismatchGapId(References & refs, int alignIdx, uint32_t & mismatches, uint32_t & gaps, uint32_t & id)
{
	const char to_char[5] = { 'A','C','G','T','N' };

	if (alignIdx >= hits_align_info.alignv.size()) return; // index exceeds the size of the alignment vector

	mismatches = 0, gaps = 0, id = 0;

	int32_t qb = hits_align_info.alignv[alignIdx].ref_begin1; // index of the first char in the reference matched part
	int32_t pb = hits_align_info.alignv[alignIdx].read_begin1; // index of the first char in the read matched part

	std::string refseq = refs.buffer[hits_align_info.alignv[alignIdx].ref_seq].sequence;

	for (uint32_t c2 = 0; c2 < hits_align_info.alignv[alignIdx].cigar.size(); ++c2)
	{
		uint32_t letter = 0xf & hits_align_info.alignv[alignIdx].cigar[c2]; // 4 low bits
		uint32_t length = (0xfffffff0 & hits_align_info.alignv[alignIdx].cigar[c2]) >> 4; // high 28 bits i.e. 32-4=28
		if (letter == 0)
		{
			for (uint32_t u = 0; u < length; ++u)
			{
				//if ((char)to_char[(int)refseq[qb]] != (char)to_char[(int)isequence[pb]]) ++mismatches;
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