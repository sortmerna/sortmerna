#pragma once
/* post-alignment tasks like calculating statistics */

#include <string>
#include <vector>

// forward
struct Runopts;
class Refstats;
struct Readstats;
class Readfeed;
class References;
class KeyValueDatabase;

void writeSummary(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts);

/**
 * Summary report (log) data structure
 */
class Summary {
public:
	// vars
	bool is_de_novo;
	bool is_otumapout;
	std::string cmd;
	std::string timestamp;
	std::string pid_str;
	uint64_t total_reads;
	uint64_t total_mapped; // total passing SW
	uint64_t total_denovo; // total SW + !ID + !COV
	uint64_t total_sw_id_cov; // total SW + ID + COV
	uint64_t total_otu; // total number of OTU groups
	uint32_t min_read_len;
	uint32_t max_read_len;
	uint64_t all_reads_len;
	std::vector<std::pair<std::string, float>> db_matches;

	// methods
	Summary();
	void write(Refstats& refstats, Readstats& readstats, Runopts& opts);
	std::string to_string(Refstats& refstats, Runopts& opts);
};
