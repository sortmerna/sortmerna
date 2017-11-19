#pragma once
/**
* FILE: alignment2.hpp
* Created: Nov 11, 2017 Mon
* Wrapper of a Reads' record and its Match results
*/

#include "read.hpp"
#include "index.hpp"
#include "references.hpp"
#include "output.hpp"

void compute_lis_alignment2(
	Read & read, Index & index, References & refs, Readstats & readstats, Output & output,
	bool & search,
	uint32_t max_SW_score,
	bool& read_to_count
);
