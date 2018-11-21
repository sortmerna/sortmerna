#pragma once
/**
 * FILE: ssw.hpp
 * Created: Oct 28, 2017 Sat
 */

#include <stdint.h>
#include <vector>
#include <iterator>

typedef struct s_align2 {
	std::vector<uint32_t> cigar;
	uint32_t ref_seq;
	int32_t ref_begin1;
	int32_t ref_end1;
	int32_t	read_begin1;
	int32_t read_end1;
	uint32_t readlen;
	uint16_t score1;
	uint16_t part;
	uint16_t index_num;
	bool strand; // flags whether this alignment was done on a forward (true) or reverse-complement (false) read.

	// default construct
	s_align2() : ref_seq(0), ref_begin1(0), ref_end1(0), read_begin1(0), read_end1(0), readlen(0), score1(0), part(0), index_num(0) {}

	// construct from binary string
	s_align2(std::string bstr)
	{
		size_t offset = 0;
		size_t len = 0; // cigar, refname length

		// cannot use copy/copy_n - VC error C4996: 'std::copy_n::_Unchecked_iterators
		std::memcpy(static_cast<void*>(&len), bstr.data(), sizeof(len));
		offset += sizeof(len);

		cigar.assign(len, 0);
		std::memcpy(static_cast<void*>(cigar.data()), bstr.data() + offset, len * sizeof(uint32_t));
		offset += len * sizeof(uint32_t);
		// ref_seq
		std::memcpy(static_cast<void*>(&ref_seq), bstr.data() + offset, sizeof(ref_seq));
		offset += sizeof(ref_seq);
		// ref_begin1
		std::memcpy(static_cast<void*>(&ref_begin1), bstr.data() + offset, sizeof(ref_begin1));
		offset += sizeof(ref_begin1);
		// ref_end1
		std::memcpy(static_cast<void*>(&ref_end1), bstr.data() + offset, sizeof(ref_end1));
		offset += sizeof(ref_end1);
		// read_begin1
		std::memcpy(static_cast<void*>(&read_begin1), bstr.data() + offset, sizeof(read_begin1));
		offset += sizeof(read_begin1);
		// read_end1
		std::memcpy(static_cast<void*>(&read_end1), bstr.data() + offset, sizeof(read_end1));
		offset += sizeof(read_end1);
		// readlen
		std::memcpy(static_cast<void*>(&readlen), bstr.data() + offset, sizeof(readlen));
		offset += sizeof(readlen);
		// score1
		std::memcpy(static_cast<void*>(&score1), bstr.data() + offset, sizeof(score1));
		offset += sizeof(score1);
		// part
		std::memcpy(static_cast<void*>(&part), bstr.data() + offset, sizeof(part));
		offset += sizeof(part);
		// index_num
		std::memcpy(static_cast<void*>(&index_num), bstr.data() + offset, sizeof(index_num));
		offset += sizeof(index_num);
		// strand
		std::memcpy(static_cast<void*>(&strand), bstr.data() + offset, sizeof(strand));
		offset += sizeof(strand);
	}

	// convert to binary string
	std::string toString()
	{
		std::string buf;
		// length of cigar
		size_t cigarlen = cigar.size();
		char* beginIt = static_cast<char*>(static_cast<void*>(&cigarlen));
		std::copy_n(beginIt, sizeof(cigarlen), std::back_inserter(buf));
		// cigar
		beginIt = static_cast<char*>(static_cast<void*>(cigar.data()));
		std::copy_n(beginIt, cigarlen * sizeof(cigar[0]), std::back_inserter(buf));
		
		beginIt = static_cast<char*>(static_cast<void*>(&ref_seq));
		std::copy_n(beginIt, sizeof(ref_seq), std::back_inserter(buf)); // ref_seq
		beginIt = static_cast<char*>(static_cast<void*>(&ref_begin1));
		std::copy_n(beginIt, sizeof(ref_begin1), std::back_inserter(buf)); // ref_begin1
		beginIt = static_cast<char*>(static_cast<void*>(&ref_end1));
		std::copy_n(beginIt, sizeof(ref_end1), std::back_inserter(buf)); // ref_end1
		beginIt = static_cast<char*>(static_cast<void*>(&read_begin1));
		std::copy_n(beginIt, sizeof(read_begin1), std::back_inserter(buf)); // read_begin1
		beginIt = static_cast<char*>(static_cast<void*>(&read_end1));
		std::copy_n(beginIt, sizeof(read_end1), std::back_inserter(buf)); // read_end1
		beginIt = static_cast<char*>(static_cast<void*>(&readlen));
		std::copy_n(beginIt, sizeof(readlen), std::back_inserter(buf)); // readlen
		beginIt = static_cast<char*>(static_cast<void*>(&score1));
		std::copy_n(beginIt, sizeof(score1), std::back_inserter(buf)); // score1
		beginIt = static_cast<char*>(static_cast<void*>(&part));
		std::copy_n(beginIt, sizeof(part), std::back_inserter(buf)); // part
		beginIt = static_cast<char*>(static_cast<void*>(&index_num));
		std::copy_n(beginIt, sizeof(index_num), std::back_inserter(buf)); // index_num
		// strand
		beginIt = static_cast<char*>(static_cast<void*>(&strand));
		std::copy_n(beginIt, sizeof(strand), std::back_inserter(buf));

		return buf;
	} // ~toString

	// for serialization
	size_t size() {
		return sizeof(uint32_t) * cigar.size()
			+ sizeof(ref_seq)
			+ sizeof(ref_begin1)
			+ sizeof(ref_end1)
			+ sizeof(read_begin1)
			+ sizeof(read_end1)
			+ sizeof(readlen)
			+ sizeof(score1)
			+ sizeof(part)
			+ sizeof(index_num)
			+ sizeof(strand);
	}

	bool operator==(const s_align2& other)
	{
		return other.cigar == cigar &&
			other.ref_seq == ref_seq &&
			other.ref_begin1 == ref_begin1 &&
			other.ref_end1 == ref_end1 &&
			other.read_begin1 == read_begin1 &&
			other.read_end1 == read_end1 &&
			other.readlen == readlen &&
			other.score1 == score1 &&
			other.part == part &&
			other.index_num == index_num &&
			other.strand == strand;
	}
} s_align2;