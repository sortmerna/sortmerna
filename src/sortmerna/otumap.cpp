/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

#include <fstream>
#include <chrono>
#include <thread>
#include <filesystem>
#include <cmath>  // std::floor

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
	INFO_NE("merging OTU map. Map vector size: ", mapv.size());
	auto starts = std::chrono::high_resolution_clock::now();
	for (unsigned i = 1; i < mapv.size(); ++i) {
		for (auto it = mapv[i].begin(); it != mapv[i].end(); ++it) {
			auto el = mapv[0].find(it->first);
			if (el == mapv[0].end())
				mapv[0].emplace(it->first, it->second);
			else
				el->second.insert(el->second.begin(), it->second.begin(), it->second.end());
		}
		mapv[i].clear();
	}
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
	INFO_NS(" ... done in [", elapsed.count(), "] sec\n");
}

void OtuMap::write()
{
	unsigned c_reads = 0;
	unsigned c_group = 0;
	if (mapv.size() > 0) {
		std::ofstream ofs;
		ofs.open(fmap);
		if (!ofs.is_open()) {
			ERR("Failed to open: ", fmap);
			exit(1);
		}

		for (auto const& amap : mapv) {
			for (auto const& apair : amap) {
				ofs << apair.first << "\t"; // ref
				unsigned i = 0;
				for (auto const& aread : apair.second) {
					ofs << aread;
					if (i < apair.second.size() - 1)
						ofs << "\t";
					++i;
					++c_reads;
				}
				ofs << std::endl;
				++c_group;
			}
		}
		if (ofs.is_open()) ofs.close();
	}
	else {
		INFO("OTU map is empty - nothing to write ");
	}
	INFO("OTU map done. Num groups: ", c_group, " num reads: ", c_reads);
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
	INFO("using OTU map file: ", fmap.generic_string());
}

size_t OtuMap::count_otu()
{
	size_t num_otu = 0;
	for (auto const& amap : mapv) {
		num_otu += amap.size();
	}
	return num_otu;
}

/*
  runs in a thread
*/
void fill_otu_map2(int id, OtuMap& otumap, Readfeed& readfeed, References& refs, KeyValueDatabase& kvdb, Runopts& opts)
{
	unsigned c_reads = 0;  // all reads count
	unsigned c_aligned = 0; // aligned reads
	unsigned c_yid_ycov = 0; // passing ID and COV
	//unsigned c_yid_ncov = 0;
	//unsigned c_nid_ycov = 0;
	//unsigned c_nid_ncov = 0;
	std::string readstr;

	if (opts.dbg_level == 2)
		INFO("OTU map thread ", id, " : ", std::this_thread::get_id(), " started");

	for (;readfeed.next(id, readstr);)
	{
		{
			Read read(readstr);
			read.init(opts);
			read.load_db(kvdb);

			if (!read.isValid)
				continue;

			if (read.is03) read.flip34();
			if (read.c_yid_ycov > 0) {
				for (auto const& align: read.alignment.alignv) {
					// process alignments that match currently loaded reference part
					if (align.index_num == refs.num && align.part == refs.part) {
						auto miss_gap_match = read.calc_miss_gap_match(refs, align);
						auto idr = std::floor(std::get<3>(miss_gap_match) * 1000.0 + 0.5) * 0.001; // round to 3 decimal
						auto covr = std::floor(std::get<4>(miss_gap_match) * 1000.0 + 0.5) * 0.001;
						auto is_id = idr >= opts.min_id;
						auto is_cov = covr >= opts.min_cov;
						if (is_id && is_cov) {
							// get reference sequence identifier
							auto refhead = refs.buffer[align.ref_num].header;
							auto ref_seq_str = refhead.substr(0, refhead.find(' '));
							// left trim '>' or '@'
							ref_seq_str.erase(ref_seq_str.begin(),
								std::find_if(ref_seq_str.begin(), ref_seq_str.end(),
									[](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));

							// get read identifier
							std::string read_seq_str = read.getSeqId();
							otumap.push(id, ref_seq_str, read_seq_str); // thread safe
							++c_yid_ycov;
						}
					}
				}
			} // ~for all alignments of a read

			readstr.resize(0);
			++c_reads;
			if (read.is_hit) ++c_aligned;
		} // ~ a read scope ends
	} // ~for all reads

	INFO("OTU map thread ", id, " : ", std::this_thread::get_id(),
		" done. All reads: ", c_reads, ". aligned reads: ", c_aligned, " c_yid_ycov: ", c_yid_ycov);
} // ~fill_otu_map2

void fill_otu_map(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts, bool is_write)
{
	INFO("==== OTU groups processing started ====");
	auto ss = std::chrono::high_resolution_clock::now();
	if (readstats.n_yid_ycov.load(std::memory_order_relaxed) > 0) {
		int numThreads = 0;
		if (opts.feed_type == FEED_TYPE::LOCKLESS)
		{
			numThreads = opts.num_read_thread_pp + opts.num_proc_thread_pp;
			INFO("using total threads: ", numThreads, 
				" including Read threads: ", opts.num_read_thread_pp, 
				" Processor threads: ", opts.num_proc_thread_pp);
		}
		else {
			numThreads = opts.num_proc_thread;
			INFO("using total threads: ", numThreads);
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
				INFO_NE("loading reference ", idx, " part ", ipart + 1, "/", refstats.num_index_parts[idx]);
				auto starts = std::chrono::high_resolution_clock::now();
				refs.load(idx, ipart, opts, refstats);
				std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
				INFO_NS(" ... done in ", elapsed.count(), " sec\n");

				starts = std::chrono::high_resolution_clock::now(); // index processing starts

				// add Readfeed job if necessary
				if (opts.feed_type == FEED_TYPE::LOCKLESS)
				{
					//tpool.addJob(f_readfeed_run);
					//tpool.addJob(Readfeed(opts.feed_type, opts.readfiles, opts.is_gz));
				}
				else if (opts.feed_type == FEED_TYPE::SPLIT_READS) {
					for (int i = 0; i < numThreads; ++i) {
						tpool.emplace_back(std::thread(fill_otu_map2, i, std::ref(otumap), 
							std::ref(readfeed), std::ref(refs),	std::ref(kvdb), std::ref(opts)));
					}
				}

				// wait till processing is done on one index part
				//tpool.waitAll(); 
				for (auto& thr: tpool) {
					thr.join();
				}

				refs.unload();
				//read_queue.reset();

				elapsed = std::chrono::high_resolution_clock::now() - starts;
				INFO_MEM("done reference ", idx, " Part: ", ipart + 1, " in ", elapsed.count(), " sec");

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
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - ss;
	INFO("==== OTU groups processing done in ", elapsed.count(), " sec ====\n");
} // ~fill_otu_map