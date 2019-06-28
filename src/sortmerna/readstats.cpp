/**
* FILE: readstats.hpp
* Created: Nov 17, 2017 Fri
* @copyright 2016-19 Clarity Genomics BVBA
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
#include <cstring> // memcpy
#include <ios>
#include <filesystem>

// 3rd party
#include "zlib.h"
#include "kseq_load.hpp"

// SMR
#include "readstats.hpp"
#include "kvdb.hpp"
#include "gzip.hpp"

// forward
std::string string_hash(const std::string &val);

Readstats::Readstats(Runopts &opts, KeyValueDatabase &kvdb)
	:
	min_read_len(MAX_READ_LEN),
	max_read_len(0),
	total_reads_mapped(0),
	total_reads_mapped_cov(0),
	all_reads_count(0),
	all_reads_len(0),
	reads_matched_per_db(opts.indexfiles.size(), 0),
	total_reads_denovo_clustering(0),
	is_stats_calc(false),
	is_total_reads_mapped_cov(false)
{
	// calculate this->dbkey
	std::string key_str_tmp("");
	for (auto readsfile : opts.readfiles)
	{
		if (key_str_tmp.size() == 0)
			key_str_tmp += std::filesystem::path(readsfile).filename().string();
		else
			key_str_tmp += "_" + std::filesystem::path(readsfile).filename().string();
	}
	dbkey = string_hash(key_str_tmp);

	bool is_restored = restoreFromDb(kvdb);

	calcSuffix(opts);

	if (!opts.exit_early)
	{
		if (!is_restored || !(is_restored && all_reads_count > 0 && all_reads_len > 0))
		{
			calculate(opts);
			store_to_db(kvdb);
		}
		else
		{
			std::cout << STAMP << "Found reads statistics in the KVDB: all_reads_count= " << all_reads_count 
				<< " all_reads_len= " << all_reads_len << " Skipping calculation..." << std::endl;
		}
	}
} // ~Readstats::Readstats

/** 
 * Go through the reads file, collect, and store in the DB the following statistics:
 *   - Total number of reads
 *   - Total length of all sequences
 *
 * Original code used a single file for paired reads => statistics is to be collected from both separate files
 */
void Readstats::calculate(Runopts &opts)
{
	std::stringstream ss;

	for (auto readfile : opts.readfiles)
	{
		std::ifstream ifs(readfile, std::ios_base::in | std::ios_base::binary);
		if (!ifs.is_open()) {
			ss << STAMP << "Failed to open Reads file: " << readfile;
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
		else
		{
			std::string line; // line from the Reads file
			std::string sequence; // full sequence of a Read (can contain multiple lines for Fasta files)
			bool isFastq = false;
			bool isFasta = false;
			uint64_t tcount = 0; // count of lines in a file
			Gzip gzip(opts.is_gz);

			auto t = std::chrono::high_resolution_clock::now();

			std::cout << STAMP << "Starting statistics calculation on file: '" << readfile << "'  ...   ";

			for (int count = 0, stat = 0; ; ++count) // count of lines in a Single record
			{
				stat = gzip.getline(ifs, line);
				++tcount;

				if (stat == RL_END)
				{
					if (!sequence.empty())
					{
						// process the last record
						++all_reads_count;
						all_reads_len += sequence.length();

						// update the minimum sequence length
						if (sequence.size() < min_read_len.load())
							min_read_len = static_cast<uint32_t>(sequence.size());

						// update the maximum sequence length
						if (sequence.size() > max_read_len.load())
							max_read_len = static_cast<uint32_t>(sequence.size());
					}
					break;
				}

				if (stat == RL_ERR)
				{
					ss.str("");
					ss << STAMP << "Failed reading from file '" << readfile << "' Exiting...";
					ERR(ss.str());
					exit(EXIT_FAILURE);
				}

				if (line.empty())
				{
					--count;
					--tcount;
					continue; // skip empty line
				}

				// right-trim whitespace in place (removes '\r' too)
				line.erase(std::find_if(line.rbegin(), line.rend(), [l = std::locale{}](auto ch) { return !std::isspace(ch, l); }).base(), line.end());
				//line.erase(std::remove_if(begin(line), end(line), [l = std::locale{}](auto ch) { return std::isspace(ch, l); }), end(line));

				if (tcount == 1)
				{
					isFastq = (line[0] == FASTQ_HEADER_START);
					isFasta = (line[0] == FASTA_HEADER_START);

					if (!(isFasta || isFastq))
					{
						ss.str("");
						ss << STAMP << " the line [" << line << "] is not FASTA/Q header";
						ERR(ss.str());
						exit(EXIT_FAILURE);
					}
				}

				if (count == 4 && isFastq)
				{
					count = 0;
					if (line[0] != FASTQ_HEADER_START)
					{
						ss.str("");
						ss << STAMP << "the line [" << line << "] is not FASTQ header. all_reads_count = "
							<< all_reads_count << " tcount= " << tcount;
						ERR(ss.str());
						exit(EXIT_FAILURE);
					}
				}

				// fastq: 0(header), 1(seq), 2(+), 3(quality)
				// fasta: 0(header), 1(seq)
				if ((isFasta && line[0] == FASTA_HEADER_START) || (count == 0 && isFastq))
				{
					if (!sequence.empty())
					{ // process previous sequence
						++all_reads_count;
						all_reads_len += sequence.length();

						// update the minimum sequence length
						if (sequence.size() < min_read_len.load())
							min_read_len = static_cast<uint32_t>(sequence.size());

						// update the maximum sequence length
						if (sequence.size() > max_read_len.load())
							max_read_len = static_cast<uint32_t>(sequence.size());
					}

					count = 0; // FASTA record start
					sequence.clear(); // clear container for the new record
				} // ~if header line
				else
				{
					if (isFastq)
					{
						if (count > 3)
						{
							ss.str("");
							ss << STAMP << " Unexpected number of lines : " << count
								<< " in a single FASTQ Read. Total reads processed: " << all_reads_count
								<< " Last sequence: " << sequence
								<< " Last line read: " << line
								<< " Exiting...";
							ERR(ss.str());
							exit(EXIT_FAILURE);
						}
						if (count == 3 || line[0] == '+')
							continue; // fastq.quality
					} // ~if fastq

					sequence += line; // fasta multiline sequence
				}
			} // ~for getline

			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
			ss.str("");
			ss << std::setprecision(2) << std::fixed << STAMP
				<< "Done statistics on file. Elapsed time: " << elapsed.count()
				<< " sec. all_reads_count= " << all_reads_count << std::endl;
			std::cout << ss.str();

			ifs.close();
			gzip.~Gzip();
		}
	} // ~for iterating reads files
} // ~Readstats::calculate

// determine the suffix (fasta, fastq, ...) of aligned strings
void Readstats::calcSuffix(Runopts &opts)
{
	const std::string suff = opts.readfiles[0].substr(opts.readfiles[0].rfind('.') + 1);
	if (suff.length() > 0 && !opts.is_gz)
		suffix.assign(suff);
	else if (filesig == FASTA_HEADER_START)
		suffix.assign("fasta");
	else
		suffix.assign("fastq");
}

/**
 * put readstats data into binary string for storing in DB
 */
std::string Readstats::toBstring()
{
	std::string buf;
	// min_read_len
	std::copy_n(static_cast<char*>(static_cast<void*>(&min_read_len)), sizeof(min_read_len), std::back_inserter(buf));
	// max_read_len
	std::copy_n(static_cast<char*>(static_cast<void*>(&max_read_len)), sizeof(max_read_len), std::back_inserter(buf));
	// total_reads_mapped (atomic int)
	auto val = total_reads_mapped.load();
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// total_reads_mapped_cov (atomic int)
	val = total_reads_mapped_cov.load();
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// all_reads_count (int)
	std::copy_n(static_cast<char*>(static_cast<void*>(&all_reads_count)), sizeof(all_reads_count), std::back_inserter(buf));
	// all_reads_len (int)
	std::copy_n(static_cast<char*>(static_cast<void*>(&all_reads_len)), sizeof(all_reads_len), std::back_inserter(buf));
	// total_reads_denovo_clustering (int)
	std::copy_n(static_cast<char*>(static_cast<void*>(&total_reads_denovo_clustering)), sizeof(total_reads_denovo_clustering), std::back_inserter(buf));
	// reads_matched_per_db (vector)
	size_t reads_matched_per_db_size = reads_matched_per_db.size();
	std::copy_n(static_cast<char*>(static_cast<void*>(&reads_matched_per_db_size)), sizeof(reads_matched_per_db_size), std::back_inserter(buf));
	for (auto entry: reads_matched_per_db)
		std::copy_n(static_cast<char*>(static_cast<void*>(&entry)), sizeof(entry), std::back_inserter(buf));

	// stats_calc_done (bool)
	std::copy_n(static_cast<char*>(static_cast<void*>(&is_stats_calc)), sizeof(is_stats_calc), std::back_inserter(buf));
	// is_total_reads_mapped_cov
	std::copy_n(static_cast<char*>(static_cast<void*>(&is_total_reads_mapped_cov)), sizeof(is_total_reads_mapped_cov), std::back_inserter(buf));

	return buf;
} // ~Readstats::toBstring

/**
 * generate human readable data representation of this object 
 */
std::string Readstats::toString()
{
	std::stringstream ss;
	ss << "min_read_len= " << min_read_len 
		<< " max_read_len= " << max_read_len
		<< " all_reads_count= " << all_reads_count
		<< " all_reads_len= " << all_reads_len
		<< " total_reads_mapped= " << total_reads_mapped
		<< " total_reads_mapped_cov= " << total_reads_mapped_cov
		<< " reads_matched_per_db= " << "TODO"
		<< " is_total_reads_mapped_cov= " << is_total_reads_mapped_cov
		<< " is_stats_calc= " << is_stats_calc << std::endl;
	return ss.str();
} // ~Readstats::toString

void Readstats::set_is_total_reads_mapped_cov()
{
	if (!is_total_reads_mapped_cov && total_reads_mapped_cov > 0)
		is_total_reads_mapped_cov = true;
}

/**
 * restore Readstats object using values stored in Key-value database 
 */
bool Readstats::restoreFromDb(KeyValueDatabase & kvdb)
{
	bool ret = false;
	size_t offset = 0;
	std::stringstream ss;

	std::string bstr = kvdb.get(dbkey);
	if (bstr.size() > 0) 
	{
		// min_read_len
		std::memcpy(static_cast<void*>(&min_read_len), bstr.data() + offset, sizeof(min_read_len));
		offset += sizeof(min_read_len);
		// max_read_len
		std::memcpy(static_cast<void*>(&max_read_len), bstr.data() + offset, sizeof(max_read_len));
		offset += sizeof(max_read_len);
		// total_reads_mapped
		uint64_t val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		total_reads_mapped = val;
		offset += sizeof(val);
		// total_reads_mapped_cov
		val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		total_reads_mapped_cov = val;
		offset += sizeof(val);
		// all_reads_count
		std::memcpy(static_cast<void*>(&all_reads_count), bstr.data() + offset, sizeof(all_reads_count));
		offset += sizeof(all_reads_count);
		// all_reads_len
		std::memcpy(static_cast<void*>(&all_reads_len), bstr.data() + offset, sizeof(all_reads_len));
		offset += sizeof(all_reads_len);
		// total_reads_denovo_clustering
		std::memcpy(static_cast<void*>(&total_reads_denovo_clustering), bstr.data() + offset, sizeof(total_reads_denovo_clustering));
		offset += sizeof(total_reads_denovo_clustering);
		// reads_matched_per_db
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
			ss.str("");
			ss << STAMP << "reads_matched_per_db.size stored in DB: " << reads_matched_per_db_size
				<< " doesn't match the number of reference files: " << reads_matched_per_db.size() << std::endl;
			std::cout << ss.str();
			ret = false;
		}

		// stats_calc_done
		std::memcpy(static_cast<void*>(&is_stats_calc), bstr.data() + offset, sizeof(is_stats_calc));
		offset += sizeof(is_stats_calc);

		// stats_calc_done
		std::memcpy(static_cast<void*>(&is_total_reads_mapped_cov), bstr.data() + offset, sizeof(is_total_reads_mapped_cov));
		offset += sizeof(is_total_reads_mapped_cov);
	} // ~if data found in DB

	return ret;
} // ~Readstats::restoreFromDb

/* push entry to Readstats::otu_map*/
void Readstats::pushOtuMap(std::string & ref_seq_str, std::string & read_seq_str)
{
	//std::lock_guard<std::mutex> omlg(otu_map_lock);
	otu_map[ref_seq_str].push_back(read_seq_str);
}

void Readstats::printOtuMap(std::string otumapfile)
{
	std::stringstream ss;
	std::ofstream omstrm;
	omstrm.open(otumapfile);

	ss << STAMP << "Printing OTU Map.." << std::endl;
	std::cout << ss.str(); ss.str("");

	for (std::map<std::string, std::vector<std::string>>::iterator it = otu_map.begin(); it != otu_map.end(); ++it)
	{
		omstrm << it->first << "\t";
		int i = 0;
		for (std::vector<std::string>::iterator itv = it->second.begin(); itv != it->second.end(); ++itv)
		{
			if (i < it->second.size() - 1)
				omstrm << *itv << "\t";
			else
				omstrm << *itv; // last element
			++i;
		}
		omstrm << std::endl;
	}
	if (omstrm.is_open()) omstrm.close();
}

void Readstats::store_to_db(KeyValueDatabase & kvdb)
{
	kvdb.put(dbkey, toBstring());

	std::cout << STAMP << "Stored Reads statistics to DB:\n    " << toString() << std::endl;
}