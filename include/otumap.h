#pragma once
/*
*/
#include <string>
#include <vector>
#include <map>

// forward
struct Runopts;
class Refstats;
struct Readstats;
class Readfeed;
class References;
class KeyValueDatabase;

class OtuMap {
	// Clustering of reads around references by similarity i.e. {ref: [read, read, ...] , ref : [read, read...] , ...}
	// calculated after alignment is done on all reads
	//  TODO: Store in DB ? Can be very big.
	std::vector<std::map<std::string, std::vector<std::string>>> mapv;
	//std::map<std::string, std::vector<std::string>> otu_map;
	//          |_Ref_ID       |_Read_IDs
public:
	std::filesystem::path fmap;
	uint64_t total_otu; // total count of OTU groups in otu_map
public:
	OtuMap(int numThreads=1);
	void push(int idx, std::string& ref_seq_str, std::string& read_seq_str);
	void merge();
	void write();
	void init(Runopts& opts);
	uint64_t count_otu();
};

void fill_otu_map(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts, bool is_write=true);
