/**
* FILE: readstats.hpp
* Created: Nov 17, 2017 Fri
* @copyright 2016-20 Clarity Genomics BVBA
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
#include "izlib.hpp"

// forward
std::string string_hash(const std::string &val); // util.cpp
std::string to_lower(std::string& val); // util.cpp

Readstats::Readstats(uint64_t all_reads_count, uint64_t all_reads_len, uint32_t min_read_len, uint32_t max_read_len, KeyValueDatabase& kvdb, Runopts& opts)
	:
	min_read_len(max_read_len),
	max_read_len(min_read_len),
	total_aligned(0),
	total_aligned_id(0),
	total_aligned_cov(0),
	total_aligned_id_cov(0),
	num_short(0),
	all_reads_count(all_reads_count),
	all_reads_len(all_reads_len),
	reads_matched_per_db(opts.indexfiles.size(), 0),
	total_denovo(0),
	total_otu(),
	is_stats_calc(false),
	is_set_aligned_id_cov(false)
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
			store_to_db(kvdb);
		}
		else
		{
			INFO("Found reads statistics in the KVDB: all_reads_count= ", all_reads_count, 
				" all_reads_len= ", all_reads_len);
		}
	}
} // ~Readstats::Readstats

// determine the suffix (fasta, fastq, ...) of aligned strings
// use the same suffix as the original reads file without 'gz' if gzipped.
void Readstats::calcSuffix(Runopts &opts)
{
	size_t pos = opts.readfiles[0].rfind('.'); // find last '.' position
	size_t pos2 = 0;
	std::string sfx = opts.readfiles[0].substr(pos + 1);
	std::string sfx_lower = to_lower(sfx);

	//if (opts.is_gz && "gz" == sfx_lower)
	if ("gz" == sfx_lower)
	{
		pos2 = opts.readfiles[0].rfind('.', pos - 1);
		sfx = opts.readfiles[0].substr(pos2 + 1, pos - pos2 - 1);
	}

	suffix.assign(sfx);
}

/**
 * put readstats data into binary string for storing in DB
 */
std::string Readstats::toBstring()
{
	std::string buf;
	// 1
	std::copy_n(static_cast<char*>(static_cast<void*>(&all_reads_count)), sizeof(all_reads_count), std::back_inserter(buf));
	// 2
	std::copy_n(static_cast<char*>(static_cast<void*>(&all_reads_len)), sizeof(all_reads_len), std::back_inserter(buf));
	// 3
	std::copy_n(static_cast<char*>(static_cast<void*>(&min_read_len)), sizeof(min_read_len), std::back_inserter(buf));
	// 4
	std::copy_n(static_cast<char*>(static_cast<void*>(&max_read_len)), sizeof(max_read_len), std::back_inserter(buf));
	// 5
	auto val = total_aligned.load(std::memory_order_relaxed);
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// 6
	val = total_aligned_id.load(std::memory_order_relaxed);
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// 7
	val = total_aligned_cov.load(std::memory_order_relaxed);
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// 8
	val = total_aligned_id_cov.load(std::memory_order_relaxed);
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// 9
	val = total_denovo.load(std::memory_order_relaxed);
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// 10
	val = num_short.load(std::memory_order_relaxed);
	std::copy_n(static_cast<char*>(static_cast<void*>(&val)), sizeof(val), std::back_inserter(buf));
	// 11
	size_t reads_matched_per_db_size = reads_matched_per_db.size();
	std::copy_n(static_cast<char*>(static_cast<void*>(&reads_matched_per_db_size)), sizeof(reads_matched_per_db_size), std::back_inserter(buf));
	// 11.1
	for (auto entry: reads_matched_per_db)
		std::copy_n(static_cast<char*>(static_cast<void*>(&entry)), sizeof(entry), std::back_inserter(buf));
	// 12
	std::copy_n(static_cast<char*>(static_cast<void*>(&is_stats_calc)), sizeof(is_stats_calc), std::back_inserter(buf));
	// 13
	std::copy_n(static_cast<char*>(static_cast<void*>(&is_set_aligned_id_cov)), sizeof(is_set_aligned_id_cov), std::back_inserter(buf));
	//
	return buf;
} // ~Readstats::toBstring

/**
 * generate human readable data representation of this object 
 */
std::string Readstats::toString()
{
	std::stringstream ss;
	ss << " all_reads_count= " << all_reads_count
		<< " all_reads_len= " << all_reads_len 
		<< " min_read_len= " << min_read_len
		<< " max_read_len= " << max_read_len
		<< " total_aligned= " << total_aligned
		<< " total_aligned_id= " << total_aligned_id
		<< " total_aligned_cov= " << total_aligned_cov
		<< " total_aligned_id_cov= " << total_aligned_id_cov
		<< " total_denovo= " << total_denovo
		<< " num_short= " << num_short
		<< " reads_matched_per_db= " << "TODO"
		<< " is_stats_calc= " << is_stats_calc
		<< " is_total_reads_mapped_cov= " << is_set_aligned_id_cov 
		<< std::endl;
	return ss.str();
} // ~Readstats::toString

void Readstats::set_is_set_aligned_id_cov()
{
	if (!is_set_aligned_id_cov && total_aligned_id_cov > 0)
		is_set_aligned_id_cov = true;
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
		// 1
		std::memcpy(static_cast<void*>(&all_reads_count), bstr.data() + offset, sizeof(all_reads_count));
		offset += sizeof(all_reads_count);
		// 2
		std::memcpy(static_cast<void*>(&all_reads_len), bstr.data() + offset, sizeof(all_reads_len));
		offset += sizeof(all_reads_len);
		// 3
		std::memcpy(static_cast<void*>(&min_read_len), bstr.data() + offset, sizeof(min_read_len));
		offset += sizeof(min_read_len);
		// 4
		std::memcpy(static_cast<void*>(&max_read_len), bstr.data() + offset, sizeof(max_read_len));
		offset += sizeof(max_read_len);
		// 5
		uint64_t val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		total_aligned = val;
		offset += sizeof(val);
		// 6
		val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		total_aligned_id = val;
		offset += sizeof(val);
		// 7
		val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		total_aligned_cov = val;
		offset += sizeof(val);
		// 8
		val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		total_aligned_id_cov = val;
		offset += sizeof(val);
		// 9
		val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		total_denovo = val;
		offset += sizeof(val);
		// 10
		val = 0;
		std::memcpy(static_cast<void*>(&val), bstr.data() + offset, sizeof(val));
		num_short = val;
		offset += sizeof(val);
		// 11
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
			WARN("reads_matched_per_db.size stored in DB: ", reads_matched_per_db_size, 
				" doesn't match the number of reference files: ", reads_matched_per_db.size());
			ret = false;
		}

		// 12
		std::memcpy(static_cast<void*>(&is_stats_calc), bstr.data() + offset, sizeof(is_stats_calc));
		offset += sizeof(is_stats_calc);

		// 13
		std::memcpy(static_cast<void*>(&is_set_aligned_id_cov), bstr.data() + offset, sizeof(is_set_aligned_id_cov));
		offset += sizeof(is_set_aligned_id_cov);
	} // ~if data found in DB

	return ret;
} // ~Readstats::restoreFromDb

void Readstats::store_to_db(KeyValueDatabase & kvdb)
{
	kvdb.put(dbkey, toBstring());
	INFO("Stored Reads statistics to DB:\n    ", toString());
}
