/**
 * FILE: ssw.hpp
 * Created: Oct 28, 2017 Sat
 */

#include <stdint.h>
#include <vector>
#include <iterator>

typedef struct s_align2 {
	uint16_t cigarLen; // need for serialization
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
	bool strand;
	s_align2() : cigarLen(0), cigar(10, 0), ref_seq(0), ref_begin1(0), ref_end1(0), read_begin1(0), read_end1(0), readlen(0), score1(0), part(0), index_num(0) {}

	// convert to binary string
	std::string toString()
	{
		size_t bufsize = sizeof(cigarLen) + cigar.size() * sizeof(uint32_t) + sizeof(ref_seq) + 4 * sizeof(ref_begin1) + sizeof(readlen) + 3 * sizeof(score1) + sizeof(strand);
		std::string buf(bufsize, 0);
		int bufidx = 0;

		// cigarlen
		char* pch = reinterpret_cast<char *>(&cigarLen);
		for (int i = 0; i < sizeof(cigarLen); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// cigar
		pch = reinterpret_cast<char *>(cigar.data());
		for (int i = 0; i < cigar.size() * sizeof(uint32_t); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// ref_seq
		pch = reinterpret_cast<char *>(&ref_seq);
		for (int i = 0; i < sizeof(ref_seq); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// ref_begin1, ref_end1, read_begin1, read_end1
		pch = reinterpret_cast<char *>(&ref_begin1);
		for (int i = 0; i < 4 * sizeof(ref_begin1); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// readlen
		pch = reinterpret_cast<char *>(&readlen);
		for (int i = 0; i < sizeof(readlen); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// score1, part, index_num
		pch = reinterpret_cast<char *>(&score1);
		for (int i = 0; i < 3 * sizeof(score1); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;
		// strand
		pch = reinterpret_cast<char *>(&strand);
		for (int i = 0; i < sizeof(strand); ++i, ++pch, ++bufidx) buf[bufidx] = *pch;

		return buf;
	} // ~toString
} s_align2; // ~typedef struct s_align2