#pragma once
#include <string>
#include <vector>

// forward
struct Runopts;
class Refstats;
struct Readstats;

/**
 * Summary report (log) data structure
 */
class Summary {
public:
	// vars
	bool is_de_novo_otu;
	bool is_otumapout;
	std::string cmd;
	std::string timestamp;
	std::string pid_str;
	uint64_t total_reads;
	uint64_t total_reads_denovo_clustering;
	uint64_t total_reads_mapped;
	uint64_t total_mapped_sw_id_cov;
	uint32_t min_read_len;
	uint32_t max_read_len;
	uint64_t all_reads_len;
	size_t total_otu;
	std::vector<std::pair<std::string, float>> db_matches;

	// methods
	Summary();
	void write(Runopts& opts, Refstats& refstats, Readstats& readstats);
	std::string to_string(Runopts& opts, Refstats& refstats);
};
