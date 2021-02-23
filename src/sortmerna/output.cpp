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

/*
* @FILE: output.cpp
* @Created: Nov 26, 2017 Sun
* @brief Object for outputting results in various formats
*/
#include "unistd.h"
#include <iomanip>
#include <fstream>
#include <cmath> // log, exp
#include <filesystem>
#include <functional> // std::ref

#include "options.hpp"
#include "output.hpp"
#include "references.hpp"
#include "readstats.hpp"
#include "processor.hpp"
#include "refstats.hpp"
#include "readsqueue.hpp"
#include "readfeed.hpp"

// forward
class Read;
struct Index;
class KeyValueDatabase;

Output::Output(Readfeed& readfeed, Runopts& opts, Readstats& readstats)
	: fastx(opts), fx_other(opts), blast(opts), denovo(opts), sam(opts), biom(opts)
{
	init(readfeed, opts, readstats);
}

//Output::~Output() {}

void Output::init(Readfeed& readfeed, Runopts& opts, Readstats& readstats)
{
	if (opts.is_fastx) fastx.init(readfeed, opts);
	if (opts.is_other) fx_other.init(readfeed, opts);
	if (opts.is_blast) blast.init(readfeed, opts);
	if (opts.is_sam) sam.init(readfeed, opts);
	if (opts.is_denovo) denovo.init(readfeed, opts);
} // ~Output::init


/*
 * called in a thread
*/
void report(int id,
	Readfeed& readfeed,
	References& refs,
	Refstats& refstats,
	KeyValueDatabase& kvdb,
	Output& output,
	Runopts& opts)
{
	size_t countReads = 0;
	size_t num_invalid = 0; // empty or invalid reads count
	//size_t denovo_n = 0; // count of denovo reads
	std::size_t num_reads = opts.is_paired ? 2 : 1;
	std::string readstr;
	std::vector<Read> reads; // two reads if paired, a single read otherwise

	INFO_MEM("Report Processor: ", id, " thread: ", std::this_thread::get_id(), " started.");
	//auto start = std::chrono::high_resolution_clock::now();

	for (bool isDone = false; !isDone;)
	{
		reads.clear();
		auto idx = id * readfeed.num_sense; // index into split_files array
		for (int i = 0; i < num_reads; ++i)
		{
			if (readfeed.next(idx, readstr))
			{
				reads.emplace_back(Read(readstr));
				reads[i].init(opts);
				reads[i].load_db(kvdb);
				readstr.resize(0);
				++countReads;
			}
			else {
				isDone = true;
			}
			if (opts.is_paired) idx ^= 1; // switch fwd-rev
		}

		if (!isDone) 
		{
			if (reads.back().isEmpty || !reads.back().isValid) {
				++num_invalid;
				continue;
			}

			// only needs one loop through all reads - reference file is not used
			if (refs.num == 0 && refs.part == 0) {
				if (opts.is_fastx)
					output.fastx.append(id, reads, opts, isDone);

				if (opts.is_other) 
					output.fx_other.append(id, reads, opts, isDone);

				if (opts.is_denovo) {
					auto is_dn = opts.is_paired ? reads[0].is_denovo || reads[1].is_denovo : reads[0].is_denovo;
					if (is_dn)
						output.denovo.append(id, reads, opts, isDone);
					//for (auto const& rd : reads) {
					//	if (!rd.is_id && !rd.is_cov) ++denovo_n;
						//if (rd.is_denovo) ++denovo_n;
					//}
				}
			}

			for (int i = 0; i < reads.size(); ++i) {
				if (opts.is_blast) output.blast.append(id, reads[i], refs, refstats, opts);
				if (opts.is_sam) output.sam.append(id, reads[i], refs, opts);
			} // ~for reads
		}
	} // ~for

	//std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start; // ~20 sec Debug/Win
	INFO_MEM("Report processor: ", id, " thread: ", std::this_thread::get_id(), " done. Processed reads: ", countReads, 
		" Invalid reads: ", num_invalid); // , " denovo count: ", denovo_n
} // ~report


// called from main. generateReports -> reportsJob
void writeReports(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	INFO("=== Report generation starts ===\n");
	auto start = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;

	int nthreads = 0;
	if (readfeed.type == FEED_TYPE::SPLIT_READS) {
		nthreads = opts.num_proc_thread;
		readfeed.init_reading(); // prepare readfeed
	}

	//ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS);
	std::vector<std::thread> tpool;
	tpool.reserve(nthreads);

	bool indb = readstats.restoreFromDb(kvdb);
	if (indb) {
		INFO("Restored Readstats from DB: ", indb);
	}

	Refstats refstats(opts, readstats);
	References refs;
	//ReadsQueue read_queue("queue_1", opts.queue_size_max, readstats.all_reads_count);
	Output output(readfeed, opts, readstats);

	if (opts.is_sam) output.sam.write_header(opts);

	// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t ref_idx = 0; ref_idx < (uint16_t)opts.indexfiles.size(); ++ref_idx)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[ref_idx]; ++idx_part)
		{
			INFO("Loading reference ", ref_idx, " part ", idx_part + 1, "/", refstats.num_index_parts[ref_idx], "  ... ");

			auto start_i = std::chrono::high_resolution_clock::now();

			refs.load(ref_idx, idx_part, opts, refstats);
			elapsed = std::chrono::high_resolution_clock::now() - start_i; // ~20 sec Debug/Win
			INFO("done in sec [", elapsed.count(), "]");

			start_i = std::chrono::high_resolution_clock::now(); // index processing starts

			// start processor
			if (opts.feed_type == FEED_TYPE::SPLIT_READS) {
				for (int i = 0; i < nthreads; ++i) {
					tpool.emplace_back(std::thread(report, i, std::ref(readfeed),
						std::ref(refs), std::ref(refstats), std::ref(kvdb), std::ref(output), std::ref(opts)));
				}
			}
			// wait till processing is done
			for (auto i = 0; i < tpool.size(); ++i) {
				tpool[i].join();
			}

			refs.unload();
			//read_queue.reset();

			elapsed = std::chrono::high_resolution_clock::now() - start_i; // index processing done
			INFO("Done reference ", ref_idx, " Part: ", idx_part + 1, " in sec: [", elapsed.count(), "].");

			start_i = std::chrono::high_resolution_clock::now();
			refs.unload();
			elapsed = std::chrono::high_resolution_clock::now() - start_i;
			INFO_MEM("References unloaded in sec: [", elapsed.count(), "].");
			tpool.clear();
			// rewind for the next index
			readfeed.rewind_in();
			readfeed.init_vzlib_in();

			if (!opts.is_blast && !opts.is_sam)	break;;
		} // ~for(idx_part)
	} // ~for(ref_idx)

	//output.closefiles();
	if (opts.is_fastx) {
		output.fastx.finish_deflate();
		output.fastx.closef();
		output.fastx.merge(readfeed.num_splits);
	}
	if (opts.is_other) {
		output.fx_other.finish_deflate();
		output.fx_other.closef();
		output.fx_other.merge(readfeed.num_splits);
	}
	if (opts.is_blast) {
		output.blast.finish_deflate();
		output.blast.closef();
		output.blast.merge(readfeed.num_splits);
	}
	if (opts.is_sam) {
		output.sam.closef();
		output.sam.merge(readfeed.num_splits);
	}
	if (opts.is_denovo) {
		output.denovo.closef();
		output.denovo.merge(readfeed.num_splits);
	}

	elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("=== Done Reports in sec [", elapsed.count(), "] ===\n");
} // ~writeReports