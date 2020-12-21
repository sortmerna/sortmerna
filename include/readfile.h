#pragma once

#include <filesystem>

struct Readfile {
	Readfile() : isFastq(false), isFasta(false), isZip(false), size(0), numreads(0) {}
	bool isFastq; // file is FASTQ
	bool isFasta; // file is FASTA
	bool isZip;   // true (compressed) | false (flat)
	unsigned numreads;  // max reads expected to be processed
	std::filesystem::path path;
	std::streampos size;
};
