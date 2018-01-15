/**
* FILE: readstats.hpp
* Created: Nov 17, 2017 Fri
*
* Read file statistics
*/

// standard
#include <chrono>
#include <algorithm> // remove_if
#include <iomanip> // output formatting
#include <locale> // isspace
#include <sstream>
#include <fstream> // std::ifstream
#include <iostream> // std::cout
#include <ios>

// 3rd party
#include "zlib.h"
#include "kseq_load.hpp"

// SMR
#include "readstats.hpp"
#include "kvdb.hpp"


void Readstats::calculate()
{
	std::stringstream ss;

	std::ifstream ifs(opts.readsfile, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open()) {
		ss << "Failed to open Reads file: " << opts.readsfile << "\n";
		std::cout << ss.str(); ss.str("");
		exit(EXIT_FAILURE);
	}
	else
	{
		std::string line;
		bool isFastq = true;

		auto t = std::chrono::high_resolution_clock::now();

		std::cout << "Readstats::calculate starting ...   ";

		for (int count = 0; std::getline(ifs, line); ) // count lines in One record
		{
			// skip empty line
			if (!line.empty())
			{
				// right-trim whitespace in place (removes '\r' too)
				line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
				//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));
				// fastq: 0(header), 1(seq), 2(+), 3(quality)
				// fasta: 0(header), 1(seq)
				if (line[0] == FASTA_HEADER_START || line[0] == FASTQ_HEADER_START)
				{
					count = 0; // record start
					isFastq = (line[0] == FASTQ_HEADER_START);
				} // ~if header line
				else {
					++count;
					if (count > 3) {
						ss << "Unexpected number of lines: " << count 
							<< " in a single Read. Total reads processed: " << number_total_read 
							<< " Exiting...\n";
						std::cout << ss.str(); ss.str("");
						exit(EXIT_FAILURE);
					}
					if (isFastq && (line[0] == '+' || count == 3)) continue; // fastq.quality

					++number_total_read;
					full_read_main += line.length();
				}
			}
		} // ~for getline

		std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
		ss << std::setprecision(2) << std::fixed 
			<< "Readstats::calculate done. Elapsed time: " << elapsed.count() 
			<< " sec. Reads processed: " << number_total_read << std::endl;
		std::cout << ss.str(); ss.str("");
	}
	ifs.close();
} // ~Readstats::calculate

bool Readstats::check_file_format()
{
	bool exit_early = false;
#ifdef HAVE_LIBZ
	// Check file format (if ZLIB supported)
	gzFile fp = gzopen(opts.readsfile.c_str(), "r");
	kseq_t *seq = kseq_init(fp);
#else
	FILE* fp = fopen(opts.readsfile.c_str(), "rb");
	kseq_t *seq = kseq_init(fileno(fp));
#endif
	int l;
	if ((l = kseq_read(seq)) >= 0)
		filesig = seq->last_char;
	else
	{
		fprintf(stderr, "  %sERROR%s: Line %d: %s unrecognized file format or empty file %s\n\n",
			startColor, endColor, __LINE__, __FILE__, opts.readsfile.c_str());
		exit_early = true;
	}
	kseq_destroy(seq);
#ifdef HAVE_LIBZ
	gzclose(fp);
#else
	fclose(fp);
#endif
	return exit_early;
} // ~Readstats::check_file_format

// determine the suffix (fasta, fastq, ...) of aligned strings
void Readstats::calcSuffix()
{
	const std::string suff = opts.readsfile.substr(opts.readsfile.rfind('.') + 1);
	if (suff.length() > 0 && !opts.have_reads_gz)
		suffix.assign(suff);
	else if (filesig == '>')
		suffix.assign("fasta");
	else
		suffix.assign("fastq");
}

std::string Readstats::toString()
{
	std::string buf;
	std::copy_n(static_cast<char*>(static_cast<void*>(&min_read_len)), sizeof(min_read_len), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&max_read_len)), sizeof(max_read_len), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&total_reads_mapped)), sizeof(total_reads_mapped), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&total_reads_mapped_cov)), sizeof(total_reads_mapped_cov), std::back_inserter(buf));
	std::copy_n(static_cast<char*>(static_cast<void*>(&total_reads_denovo_clustering)), sizeof(total_reads_denovo_clustering), std::back_inserter(buf));

	// vector reads_matched_per_db
	size_t reads_matched_per_db_size = reads_matched_per_db.size();
	std::copy_n(static_cast<char*>(static_cast<void*>(&reads_matched_per_db_size)), sizeof(reads_matched_per_db_size), std::back_inserter(buf));
	for (auto entry: reads_matched_per_db)
		std::copy_n(static_cast<char*>(static_cast<void*>(&entry)), sizeof(entry), std::back_inserter(buf));

	return buf;
}

bool Readstats::restoreFromDb(KeyValueDatabase & kvdb)
{
	bool ret = false;
	std::string bstr = kvdb.get("Readstats");
	if (bstr.size() == 0) { return ret; }
	size_t offset = 0;
	std::stringstream ss;

	std::memcpy(static_cast<void*>(&min_read_len), bstr.data() + offset, sizeof(min_read_len));
	offset += sizeof(min_read_len);

	std::memcpy(static_cast<void*>(&max_read_len), bstr.data() + offset, sizeof(max_read_len));
	offset += sizeof(max_read_len);

	std::memcpy(static_cast<void*>(&total_reads_mapped), bstr.data() + offset, sizeof(total_reads_mapped));
	offset += sizeof(total_reads_mapped);

	std::memcpy(static_cast<void*>(&total_reads_mapped_cov), bstr.data() + offset, sizeof(total_reads_mapped_cov));
	offset += sizeof(total_reads_mapped_cov);

	std::memcpy(static_cast<void*>(&total_reads_denovo_clustering), bstr.data() + offset, sizeof(total_reads_denovo_clustering));
	offset += sizeof(total_reads_denovo_clustering);

	size_t reads_matched_per_db_size = 0;
	std::memcpy(static_cast<void*>(&reads_matched_per_db_size), bstr.data() + offset, sizeof(reads_matched_per_db_size));
	offset += sizeof(reads_matched_per_db_size);
	if (reads_matched_per_db_size == reads_matched_per_db.size()) 
	{
		for (std::vector<uint64_t>::iterator it = reads_matched_per_db.begin(); it != reads_matched_per_db.end(); ++it)
		{
			std::memcpy(static_cast<void*>(&*it), bstr.data() + offset, sizeof(uint64_t));
			offset += sizeof(uint64_t);
		}
		ret = true;
	}
	else
	{
		ss << "Readstats::restoreFromDb: reads_matched_per_db.size stored in DB: " << reads_matched_per_db_size 
			<< " doesn't match the number of reference files: "	<< reads_matched_per_db.size() << std::endl;
		std::cout << ss.str(); ss.str("");
		ret = false;
	}

	return ret;
} // ~Readstats::restoreFromDb

void Readstats::pushOtuMap(std::string & ref_seq_str, std::string & read_seq_str)
{
	//std::lock_guard<std::mutex> omlg(otu_map_lock);
	otu_map[ref_seq_str].push_back(read_seq_str);
}

void Readstats::increment_total_reads_mapped_cov()
{
	//std::lock_guard<std::mutex> rmcg(total_reads_mapped_cov_lock);
	++total_reads_mapped_cov;
}