/*
 * FILE: references.cpp
 * Created: Dec 23, 2017 Sat
 */
#include <fstream>
#include <sstream>
#include <iostream>
#include <cctype> // std::isspace
#include <ios>
#include <cstdint>
#include <locale>

#include "references.hpp"
#include "refstats.hpp"
#include "options.hpp"
#include "common.hpp"

// prototype: 'load_ref'
void References::load(uint32_t idx_num, uint32_t idx_part, Runopts & opts, Refstats & refstats)
{
	std::stringstream ss;
	num = idx_num;
	part = idx_part;
	uint32_t numseq_part = refstats.index_parts_stats_vec[idx_num][idx_part].numseq_part;

	std::ifstream ifs(opts.indexfiles[idx_num].first, std::ios_base::in | std::ios_base::binary); // open reference file

	if (!ifs.is_open())
	{
		ss << "  " << RED << "ERROR" << COLOFF  << ": [Line " << __LINE__ << ": " << __FILE__ 
			<< "] could not open file " << opts.indexfiles[idx_num].first  << std::endl;
		std::cerr << ss.str(); ss.str("");
		exit(EXIT_FAILURE);
	}

	// set the file pointer to the first sequence added to the index for this index file section
	ifs.seekg(refstats.index_parts_stats_vec[idx_num][idx_part].start_part);
	if (ifs.fail())
	{
		ss << "  " << RED << "ERROR" << COLOFF << ": [Line " << __LINE__ << ": " << __FILE__
			<< "] could not locate the sequences used to construct the index" << std::endl
			<< "  Check that your --ref <FASTA file, index name> correspond correctly for the FASTA file: " 
			<< opts.indexfiles[idx_num].first << std::endl;
		std::cerr << ss.str(); ss.str("");
	}

	// load references sequences, skipping the empty lines & spaces
	uint64_t num_seq_read = 0;
	std::string line;
	References::BaseRecord rec;
	bool isFastq = true;
	bool lastRec = false;

	for (int count = 0; num_seq_read != numseq_part; )
	{
		if (!lastRec) std::getline(ifs, line);

		if (line.empty() && !lastRec)
		{
			if (ifs.eof()) lastRec = true;
			continue;
		}

		if (lastRec)
		{
			if (!rec.isEmpty)
			{
				buffer.push_back(rec);
				num_seq_read++;
			}
			break;
		}

		// remove whitespace (removes '\r' too)
		line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
		//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));
		//if (line[line.size() - 1] == '\r') line.erase(line.size() - 1); // remove trailing '\r'
		// fastq: 0(header), 1(seq), 2(+), 3(quality)
		// fasta: 0(header), 1(seq)
		if (line[0] == FASTA_HEADER_START || line[0] == FASTQ_HEADER_START)
		{
			if (!rec.isEmpty)
			{
				buffer.push_back(rec); // push record created before this current header
				rec.isEmpty = true;
				num_seq_read++;
				count = 0;
			}

			// start new record
			rec.clear();
			isFastq = (line[0] == FASTQ_HEADER_START);
			rec.format = isFastq ? Format::FASTQ : Format::FASTA;
			rec.header = line;
			rec.isEmpty = false;
		} // ~header or last record
		else
		{
			if (isFastq && count > 3) 
			{
				ss << "  " << RED << "ERROR" << COLOFF << " too many lines (> 4) for FASTQ file" << std::endl;
				std::cerr << ss.str(); ss.str("");
				exit(EXIT_FAILURE);
			}

			++count; // count the four FASTQ lines

			if (isFastq && line[0] == '+') continue;

			if (isFastq && count == 3)
			{
				rec.quality = line;
				continue;
			}

			convert_fix(line);
			rec.sequence += line;
		} // ~not header
		if (ifs.eof()) lastRec = true; // push and break
	} // ~for
} // ~References::load

  // convert sequence to numerical form and fix ambiguous chars
void References::convert_fix(std::string & seq)
{
	for (std::string::iterator it = seq.begin(); it != seq.end(); ++it)
	{
		if (*it != 32) // space
			*it = nt_table[(int)*it];
	}
}

void References::clear()
{
	buffer.clear(); // TODO: is this enough?
} // ~References::clear
