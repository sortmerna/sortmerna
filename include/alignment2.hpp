#pragma once
/**
* FILE: alignment2.hpp
* Created: Nov 11, 2017 Mon
* Wrapper of a Reads' record and its Match results
*/

// forward
class Read;
struct Runopts;
struct Index;
class References;
class Output;
struct Readstats;
class Refstats;

void compute_lis_alignment(
	Read & read, Runopts & opts, Index & index, References & refs, Readstats & readstats, Refstats & refstats, Output & output,
	bool & search,
	uint32_t max_SW_score,
	bool& read_to_count
);
