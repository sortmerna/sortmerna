#pragma once

#include <string>

class BaseRecord
{
public:
	std::string header;
	std::string sequence;
	std::string quality; // "" (fasta) | "xxx..." (fastq)
	std::string format = "fasta"; // fasta | fastq

	BaseRecord() {}
	~BaseRecord() {}
};
