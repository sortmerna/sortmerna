/* 
 * FILE: cmd.cpp
 * Created: Jan 17, 2018 Wed
 * @copyright 2016-20 Clarity Genomics BVBA
 */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cctype> // isdigit
#include <string>
#include <utility> // std::pair

#include "cmd.hpp"
#include "options.hpp"
#include "kvdb.hpp"
#include "read.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "references.hpp"
#include "index.hpp"
#include "reader.hpp"

const std::string OPT_DB   = "--db";
//const std::string OPT_IDX  = "--idx";
const std::string OPT_PART = "--part";
const std::string OPT_READ = "--read";
const std::string OPT_POS  = "--pos";

void CmdSession::run(Runopts & opts)
{
	std::stringstream ss;
	std::string cmd;

	for (;;)
	{
		std::cout << "Enter command: [read --id=ID | --db, index --idx=N --part=N --read=N --pos=0 | exit]: ";
		std::getline(std::cin, cmd);
		std::cout << "Processing command: " << cmd << std::endl;
		std::istringstream iss(cmd);
		std::vector<std::string> cmdv((std::istream_iterator<std::string>(iss)),
			std::istream_iterator<std::string>());
		if ("exit"  == cmdv[0]) break;
		if ("read"  == cmdv[0]) cmdRead(opts, cmd);
		if ("index" == cmdv[0]) cmdIndex(opts, cmd);
		if ("ref"   == cmdv[0]) cmd_max_ref_part(opts, cmd);
		if ("test"  == cmdv[0]) cmdTest(opts, cmd);
	}
} // ~CmdSession::run

bool getOpt(std::string & cmd, const std::string OPT, std::string & optval)
{
	std::string::size_type pos;
	bool isopt = false;

	pos = cmd.find(OPT); // option start
	if (pos != std::string::npos)
	{
		std::string::size_type pos2 = cmd.find(" ", pos); // option end
		if (pos2 == std::string::npos)
			optval = cmd.substr(pos + OPT.size() + 1);
		else
			optval = cmd.substr(pos + OPT.size() + 1, pos2 - pos - OPT.size() - 1);

		isopt = true;
	}
	else
	{
		std::cout << "getOpt: missing " << OPT << " option" << std::endl;
		isopt = false;
	}

	return isopt;
} // ~getOptval


/* 
 * @param cmdv[1]  read index e.g. 0, 10, 480 i.e. read number in the file
 */
void CmdSession::cmdRead(Runopts & opts, std::string & cmd)
{
	std::stringstream ss;

	std::string readid;
	bool isnumber = false;
	bool isdb = false; // --db load read from DB
	Read read;

	getOpt(cmd, OPT_ID, readid);
	isdb = cmd.find(OPT_DB) != std::string::npos;

	isnumber = !readid.empty() && std::find_if(readid.begin(), readid.end(),
		[](char c) { return !std::isdigit(c); }) == readid.end();

	if ( isnumber )
	{
		if (isdb)
		{
			KeyValueDatabase kvdb(opts.kvdbdir.string());
			read.clear();
			read.init(opts); // TODO: pass the required reads file number i.e. 0 or 1 to generate a correct read.id
			ss << read.matchesToJson() << std::endl;
		}
		else
		{
			read.id = readid;
			bool isok = Reader::loadReadByIdx(read);
			ss << "Read load OK " << isok << std::endl;
		}
	}
	std::cout << ss.str(); ss.str("");
} // ~CmdSession::cmdRead

/* 
 * index --idx=0 --part=1 --ref=2095 --read=480 --pos=0
 *
 * Find if the kmer at position '--pos' in the read '--read' has matches in the reference '--ref'
 */
void CmdSession::cmdIndex(Runopts & opts, std::string & cmd)
{
	std::stringstream ss;
	std::string idxval;
	std::string partval;
	std::string readid;
	std::string posval;
	std::string refid;
	bool isok = true;

	isok = isok && getOpt(cmd, OPT_IDX, idxval);
	isok = isok && getOpt(cmd, OPT_PART, partval);
	isok = isok && getOpt(cmd, OPT_READ, readid);
	isok = isok && getOpt(cmd, OPT_POS, posval);
	isok = isok && getOpt(cmd, OPT_REF, refid);

	if (!isok)
	{
		std::cout << "cmdIndex: missing some options. Returning.." << std::endl;
		return;
	}

	KeyValueDatabase kvdb(opts.kvdbdir.string());
	Readstats readstats(opts, kvdb);
	Refstats refstats(opts, readstats);
	References refs;
	Index index(opts);
	Read read;

	// find half-kmer prefix/suffix matches
	//
	read.id = readid;
	isok = Reader::loadReadByIdx(read);
	if (read.sequence.size() > 0 && read.isequence.size() == 0)
		read.seqToIntStr();
	index.load(std::stoi(idxval), std::stoi(partval), opts, refstats);
	refs.load(std::stoi(idxval), std::stoi(partval), opts, refstats);
	// find kmer prefix hash
	uint32_t kmerhash = read.hashKmer(std::stoi(posval), 9);
	if (kmerhash > index.lookup_tbl.size() - 1)
	{
		std::cout << "Hash: " << kmerhash << " is larger than Lookup table size: " << index.lookup_tbl.size() << std::endl;
		return;
	}
	std::cout << "read.id: " << readid << " Kmer position: " << posval << " DB matches: " << index.lookup_tbl[kmerhash].count << std::endl;
	//
	// find full kmer matches - burst-trie search
	//
	std::vector<UCHAR> bitvec; // window (prefix/suffix) bitvector
	uint32_t bitvec_size = (refstats.partialwin[index.index_num] - 2) << 2; // e.g. 9 - 2 = 0000 0111 << 2 = 0001 1100 = 28
	uint32_t offset = (refstats.partialwin[index.index_num] - 3) << 2; // e.g. 9 - 3 = 0000 0110 << 2 = 0001 1000 = 24
	bitvec.resize(bitvec_size);
	std::fill(bitvec.begin(), bitvec.end(), 0);

	init_win_f(&read.isequence[std::stoi(posval) + refstats.partialwin[index.index_num]],
		&bitvec[0],
		&bitvec[4],
		refstats.numbvs[index.index_num]);

	bool accept_zero_kmer = false;
	std::vector<id_win> id_hits;

	// search burst-trie
	traversetrie_align(
		index.lookup_tbl[kmerhash].trie_F,
		0,
		0,
		&bitvec[0],
		&bitvec[offset],
		accept_zero_kmer,
		id_hits,
		//read.id,
		std::stoi(posval),
		refstats.partialwin[index.index_num],
		opts
	);

	// map of k-mer occurrences on the references i.e. 
	// <reference number : number of the k-mer occurrences>
	std::map<uint32_t, uint32_t> seq_kmer_freq_map;
	std::vector<std::pair<uint32_t, uint32_t>> seq_kmer_freq_vec;

	for (auto it = id_hits.begin(); it != id_hits.end(); ++it)
	{
		// sort matches by Reference ID
		std::sort(&index.positions_tbl[it->id].arr[0],
			&index.positions_tbl[it->id].arr[index.positions_tbl[it->id].size - 1],
			[](seq_pos a, seq_pos b) { return a.seq > b.seq; });

		std::cout << "kmer iD: " << it->id << " Num hits: " << index.positions_tbl[it->id].size << std::endl;

		for ( uint32_t i = 0; i < index.positions_tbl[it->id].size; ++i)
		{
			// populate frequency map
			auto map_it = seq_kmer_freq_map.find(index.positions_tbl[it->id].arr[i].seq);
			if (map_it != seq_kmer_freq_map.end())
				map_it->second++; // increment the frequency
			else
				seq_kmer_freq_map[index.positions_tbl[it->id].arr[i].seq] = 1; // add seq to map with freq = 1

			if (index.positions_tbl[it->id].arr[i].seq == std::stoi(refid))
				std::cout << "Found match in Ref: " << std::stoi(refid) 
				<< " at Ref pos: " << index.positions_tbl[it->id].arr[i].pos 
				<< " hit number: " << i << std::endl;
		}
		//std::cout << "Max Reference number: " << index.positions_tbl[it->id].arr[0].seq << std::endl;
	}

	// copy frequency map pairs to vector
	for (auto freq_pair : seq_kmer_freq_map)
		seq_kmer_freq_vec.push_back(freq_pair);

	// sort frequency vector. apparently Map cannot be sorted
	std::sort(seq_kmer_freq_vec.begin(), seq_kmer_freq_vec.end(),
		[](std::pair<uint32_t, uint32_t> e1, std::pair<uint32_t, uint32_t> e2) {
		if (e1.second == e2.second)
			return e1.first > e2.first; // order seq numbers descending if frequencies equals
		else
			return e1.second > e2.second; // order frequencies descending
	});

	auto map_it = seq_kmer_freq_map.find(std::stoi(refid));
	if (map_it != seq_kmer_freq_map.end())
		std::cout << "Read: " << readid << " at position: " << posval << " has " << map_it->second << " matches in reference: " << refid << std::endl;
	else
		std::cout << "Read: " << readid << " at position: " << posval << " has no matches in reference: " << refid << std::endl;

} // ~CmdSession::cmdIndex

/* 
 * Get Max reference number given the reference part 
 * ref --idx=0 --part=1
 */
void CmdSession::cmd_max_ref_part(Runopts & opts, std::string & cmd)
{
	std::stringstream ss;
	std::string idxval;
	std::string partval;
	bool isok = true;

	isok = isok && getOpt(cmd, OPT_IDX, idxval);
	isok = isok && getOpt(cmd, OPT_PART, partval);

	if (!isok)
	{
		std::cout << "cmdIndex: missing some options. Returning.." << std::endl;
		return;
	}

	KeyValueDatabase kvdb(opts.kvdbdir.string());
	Readstats readstats(opts, kvdb);
	Refstats refstats(opts, readstats);
	References refs;

	refs.load(std::stoi(idxval), std::stoi(partval), opts, refstats);

	std::cout << " Reference file number: " << idxval 
		<< " Reference part: " << partval
		<< " Part size: " << refs.buffer.size()
		<< " Max Ref ID: " << refs.buffer[refs.buffer.size() - 1].id
		<< " Max Ref NID: " << refs.buffer[refs.buffer.size() - 1].nid
		<< std::endl;
 
} // ~CmdSession::cmd_max_ref_part

/* Sample command: test 
*/
void CmdSession::cmdTest(Runopts & opts, std::string & cmd)
{
	std::stringstream ss;
} // ~CmdSession::cmdTest