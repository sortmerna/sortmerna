#include <fstream>
#include <chrono>
#include <thread>
#include <filesystem>

#include "common.hpp"
#include "otumap.h"
#include "read.hpp"
#include "readfeed.hpp"
#include "references.hpp"
#include "refstats.hpp"
#include "readstats.hpp"

OtuMap::OtuMap(int numThreads) : mapv(numThreads), total_otu(0) {}

void OtuMap::push(int idx, std::string& ref_seq_str, std::string& read_seq_str)
{
	mapv[idx][ref_seq_str].emplace_back(read_seq_str);
}

void OtuMap::merge()
{
	INFO("=== Merging OTU map. Map vector size: ", mapv.size(), " ===");
	auto starts = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < mapv.size(); ++i) {
		if (i > 0) {
			mapv[0].merge(mapv[i]);
			mapv[i].clear();
		}
	}
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO("=== Merging OTU map done in [", elapsed.count(), "] sec ===\n");
}

void OtuMap::write()
{
	if (mapv.size() > 0) {
		std::ofstream ofs;
		ofs.open(fmap);
		if (!ofs.is_open()) {
			ERR("Failed to open: ", fmap);
			exit(1);
		}

		INFO("Printing OTU Map ...");
		for (std::map<std::string, std::vector<std::string>>::iterator it = mapv[0].begin(); it != mapv[0].end(); ++it)
		{
			ofs << it->first << "\t";
			int i = 0;
			for (std::vector<std::string>::iterator itv = it->second.begin(); itv != it->second.end(); ++itv)
			{
				if (i < it->second.size() - 1)
					ofs << *itv << "\t";
				else
					ofs << *itv; // last element
				++i;
			}
			ofs << std::endl;
		}
		if (ofs.is_open()) ofs.close();
	}
	else {
		INFO("OTU map is empty - nothing to write ");
	}
}

void OtuMap::init(Runopts& opts)
{
	// OTU map output file  WORKDIR/out/otu_map_PID.txt
	std::ofstream ofs;
	std::filesystem::path pdir = opts.aligned_pfx.has_parent_path() ? opts.aligned_pfx.parent_path() : opts.aligned_pfx;
	std::string bname = "otu_map";
	std::string sfx = opts.is_pid ? "_" + std::to_string(getpid()) : "";
	std::string ext = ".txt";
	fmap = pdir / (bname + sfx + ext);
	INFO("Using OTU map file: ", fmap.generic_string());
}

uint64_t OtuMap::count_otu()
{
	uint64_t num_otu = 0;
	for (int i = 0; i < mapv.size(); ++i) {
		num_otu += mapv[i].size();
	}
	return num_otu;
}

/*
  runs in a thread
*/
void fill_otu_map2(int id, OtuMap& otumap, Readfeed& readfeed, References& refs, Refstats& refstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	unsigned countReads = 0;
	unsigned count_reads_aligned = 0;
	std::string readstr;

	INFO("OTU map processor: ", id, " thread: ", std::this_thread::get_id(), " started");

	for (;readfeed.next(id, readstr);)
	{
		{
			Read read(readstr);
			read.init(opts);
			read.load_db(kvdb);

			if (!read.isValid)
				continue;

			{
				if (read.is03) read.flip34();
				if (read.is_id && read.is_cov) {
					for (int i = 0; i < read.alignment.alignv.size(); ++i) {
						// process alignments that match currently loaded reference part
						auto is_right_ref = read.alignment.alignv[i].index_num == refs.num && read.alignment.alignv[i].part == refs.part;
						if (is_right_ref) {
							// reference sequence identifier for mapped read
							std::string refhead = refs.buffer[read.alignment.alignv[i].ref_num].header;
							std::string ref_seq_str = refhead.substr(0, refhead.find(' '));
							// left trim '>' or '@'
							ref_seq_str.erase(ref_seq_str.begin(),
								std::find_if(ref_seq_str.begin(), ref_seq_str.end(),
									[](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));

							// read identifier
							std::string read_seq_str = read.getSeqId();
							otumap.push(id, ref_seq_str, read_seq_str); // thread safe
						}
					}
				}
			}
			readstr.resize(0);
			++countReads;
			if (read.is_hit) ++count_reads_aligned;
		}
	}

	INFO("OTU map processor: ", id, " thread: ", std::this_thread::get_id(),
		" done. Processed reads: ", countReads, ". count_reads_aligned: ", count_reads_aligned);
} // ~fill_otu_map2

void fill_otu_map(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts, bool is_write)
{
	INFO("==== OTU groups processing started ====");
	if (readstats.total_aligned_id_cov.load(std::memory_order_relaxed) > 0) {
		int numThreads = 0;
		if (opts.feed_type == FEED_TYPE::LOCKLESS)
		{
			numThreads = opts.num_read_thread_pp + opts.num_proc_thread_pp;
			INFO("using total threads: ", numThreads, " including Read threads: ", opts.num_read_thread_pp, " Processor threads: ", opts.num_proc_thread_pp);
		}
		else {
			numThreads = opts.num_proc_thread;
			INFO("Using total threads: ", numThreads);
			readfeed.init_reading(); // prepare readfeed
		}

		std::vector<std::thread> tpool;
		tpool.reserve(numThreads);

		Refstats refstats(opts, readstats);
		OtuMap otumap(numThreads);
		References refs;
		// loop through every reference file part
		for (uint16_t idx = 0; idx < opts.indexfiles.size(); ++idx) {
			for (uint16_t ipart = 0; ipart < refstats.num_index_parts[idx]; ++ipart) {
				// load reference
				INFO("Loading reference ", idx, " part ", ipart + 1, "/", refstats.num_index_parts[idx], "  ... ");
				auto starts = std::chrono::high_resolution_clock::now();
				refs.load(idx, ipart, opts, refstats);
				std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
				INFO("done. Elapsed sec: [", elapsed.count(), "]");

				starts = std::chrono::high_resolution_clock::now(); // index processing starts

				// add Readfeed job if necessary
				if (opts.feed_type == FEED_TYPE::LOCKLESS)
				{
					//tpool.addJob(f_readfeed_run);
					//tpool.addJob(Readfeed(opts.feed_type, opts.readfiles, opts.is_gz));
				}
				else if (opts.feed_type == FEED_TYPE::SPLIT_READS) {
					for (int i = 0; i < numThreads; ++i) {
						tpool.emplace_back(std::thread(fill_otu_map2, i, std::ref(otumap), std::ref(readfeed), std::ref(refs),
							std::ref(refstats), std::ref(kvdb), std::ref(opts)));
					}
				}

				// wait till processing is done on one index part
				//tpool.waitAll(); 
				for (auto i = 0; i < tpool.size(); ++i) {
					tpool[i].join();
				}

				refs.unload();
				//read_queue.reset();

				elapsed = std::chrono::high_resolution_clock::now() - starts;
				INFO_MEM("Done reference ", idx, " Part: ", ipart + 1, " Elapsed sec: ", elapsed.count());

				refs.unload();
				INFO_MEM("References unloaded.");
				tpool.clear();
				// rewind for the next index
				readfeed.rewind_in();
				readfeed.init_vzlib_in();
			} // ~for(ipart)
		} // ~for(idx)

		otumap.merge();
		readstats.total_otu = otumap.count_otu();
		if (is_write) {
			otumap.init(opts); // prepare the file
			otumap.write();
		}
	}
	else {
		INFO("No OTU groups to output - No reads pass %ID and %COV thresholds");
	}
	INFO("==== OTU groups processing done ====\n");
} // ~fill_otu_map