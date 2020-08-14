/**
 * FILE: processor.cpp
 * Created: Nov 26, 2017 Sun
 *
 * Callable objects designed to be run in threads
 *
 * @copyright 2016-20 Clarity Genomics BVBA
 */

#include <iostream>
#include <sstream>
#include <chrono>
#include <iomanip> // std::setprecision

#include "processor.hpp"
#include "readsqueue.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "output.hpp"
#include "index.hpp"
#include "references.hpp"
#include "options.hpp"
#include "read.hpp"
#include "ThreadPool.hpp"
#include "readsfile.hpp"

// forward
void computeStats(Read & read, Readstats & readstats, Refstats & refstats, References & refs, Runopts & opts);

/* Runs in a thread. Pops reads from the Reads Queue */
void Processor::run()
{
	unsigned num_all = 0; // all reads this processor sees
	unsigned num_skipped = 0; // reads already processed i.e. results found in Database
	unsigned num_hit = 0; // count of reads with read.hit = true found by a single thread - just for logging
	std::string readstr;
	
	INFO("Processor ", id, " thread ", std::this_thread::get_id(), " started");

	while (readQueue.pop(readstr))
	{
		{
			Read read(readstr);
			read.init(opts);
			read.is_too_short = read.sequence.size() < refstats.lnwin[index.index_num];

			if (read.is_too_short) {
				read.isValid = false;
				readstats.short_reads_num.fetch_add(1, std::memory_order_relaxed);
			}

			if (read.isValid) {
				read.load_db(kvdb);
			}

			if (read.isEmpty || !read.isValid || read.is_aligned) {
				if (read.is_aligned) {
					++num_skipped;
				}
				//INFO("Skpping read ID: ", read.id);
				continue;
			}

			// search the forward and/or reverse strands depending on Run options
			auto num_strands = 0;
			//opts.forward = true; // TODO: this discards the possiblity of forward = false
			bool search_single_strand = opts.is_forward ^ opts.is_reverse; // search only a single strand
			if (search_single_strand)
				num_strands = 1; // only search the forward xor reverse strand
			else
				num_strands = 2; // search both strands. The default when neither -F or -R were specified

			//                                                  |- stop if read was aligned on FWD strand
			for (auto count = 0; count < num_strands && !read.is_aligned; ++count)
			{
				if ((search_single_strand && opts.is_reverse) || count == 1)
				{
					if (!read.reversed)
						read.revIntStr();
				}
				// 'paralleltraversal.cpp::align_cb'
				callback(opts, index, refs, readstats, refstats, read, search_single_strand || count == 1);
				read.id_win_hits.clear(); // bug 46
			}

			// write to DB - thread safe
			if (read.isValid && !read.isEmpty)
			{
				if (read.is_hit) ++num_hit;
				if (read.is_new_hit)
					kvdb.put(read.id, read.toBinString());
			}

			readstr.resize(0);
			++num_all;
		} // ~if & read destroyed
	} // ~while there are reads

	INFO("Processor ", id, " thread ", std::this_thread::get_id(), " done. Processed ", num_all, 
		" reads. Skipped already processed: ", num_skipped, " reads", " Aligned reads (passing E-value): ", num_hit);
} // ~Processor::run

void PostProcessor::run()
{
	unsigned countReads = 0;
	unsigned count_reads_aligned = 0;
	std::string readstr;

	INFO("PostProcessor: ", id, " thread: ", std::this_thread::get_id(), " started");

	for (;readQueue.pop(readstr);)
	{
		{
			Read read(readstr);
			read.init(opts);
			read.load_db(kvdb);

			if (!read.isValid)
				continue;

			callback(read, readstats, refstats, refs, opts); // callbacks.cpp::computeStats
			readstr.resize(0);
			++countReads;
			if (read.is_hit) ++count_reads_aligned;
		}
	}

	INFO("PostProcessor: ", id, " thread: ", std::this_thread::get_id(), 
		" done. Processed reads: ",	countReads, ". count_reads_aligned: ", count_reads_aligned);

} // ~PostProcessor::run

void ReportProcessor::run()
{
	unsigned countReads = 0;
	unsigned num_invalid = 0; // empty or invalid reads count
	std::size_t num_reads = opts.is_paired ? 2 : 1;
	std::string readstr;
	std::vector<Read> reads; // two reads if paired, a single read otherwise

	INFO_MEM("Report Processor: ", id, " thread: ", std::this_thread::get_id(), " started.");

	for (bool isDone = false; !isDone;)
	{
		reads.clear();
		for (std::size_t i = 0; i < num_reads; ++i)
		{
			if (readQueue.pop(readstr))
			{
				Read read(readstr);
				read.init(opts);
				read.load_db(kvdb);
				reads.push_back(read);
				readstr.resize(0);
				++countReads;
			}
			else {
				isDone = true;
			}
		}

		if (!isDone) {
			if (reads.back().isEmpty || !reads.back().isValid) {
				++num_invalid;
			}
			else {
				callback(reads, opts, refs, refstats, output);
			}
		}
	} // ~for

	INFO_MEM("Report Processor: ", id, " thread: ", std::this_thread::get_id(), " done. Processed reads: ", countReads, " Invalid reads: ", num_invalid);
} // ~ReportProcessor::run

// called from main
void postProcess(Runopts& opts, Readstats& readstats, Output& output, KeyValueDatabase& kvdb)
{
	int N_READ_THREADS = opts.num_read_thread_pp;
	int N_PROC_THREADS = opts.num_proc_thread_pp; // opts.num_proc_threads
	int loopCount = 0; // counter of total number of processing iterations. TODO: no need here?
	
	INFO("==== Starting Post-processing (alignment statistics report) ====\n\n");

	ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS + opts.num_write_thread);
	ReadsQueue read_queue("queue_1", opts.queue_size_max, readstats.all_reads_count);
	bool indb = readstats.restoreFromDb(kvdb);

	if (indb) {	INFO("Restored Readstats from DB:\n    ", readstats.toString()); }

	readstats.total_reads_denovo_clustering = 0; // TODO: to prevent incrementing the stored value. Change this if ever using 'stats_calc_done"

	//if (!readstats.stats_calc_done)
	//{
		Refstats refstats(opts, readstats);
		References refs;

		// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
		for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); ++index_num)
		{
			// iterate parts of reference files
			for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
			{
				INFO("Loading reference ", index_num, " part ", idx_part + 1, "/", refstats.num_index_parts[index_num], "  ... ");

				auto starts = std::chrono::high_resolution_clock::now(); // index loading start
				refs.load(index_num, idx_part, opts, refstats);
				std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;

				INFO("done. Elapsed sec: [", elapsed.count(), "]");

				starts = std::chrono::high_resolution_clock::now(); // index processing starts

				// start Reader
				tpool.addJob(Readsfile(opts.readfiles, opts.is_gz));

				// start Processor
				tpool.addJob(PostProcessor("postproc_1", read_queue, opts, refs, readstats, refstats, kvdb, computeStats));

				tpool.waitAll(); // wait till processing is done on one index part

				refs.unload();
				read_queue.reset();
				++loopCount;

				elapsed = std::chrono::high_resolution_clock::now() - starts;
				INFO_MEM("Done reference ", index_num, " Part: ", idx_part + 1, " Elapsed sec: ", elapsed.count());
			} // ~for(idx_part)
		} // ~for(index_num)

		INFO("total_reads_denovo_clustering = " , readstats.total_reads_denovo_clustering);

		readstats.set_is_total_mapped_sw_id_cov();
		readstats.is_stats_calc = true;
		readstats.store_to_db(kvdb); // store reads statistics computed by post-processor
	//} // ~if !readstats.stats_calc_done

	output.writeLog(opts, refstats, readstats);

	if (opts.is_otu_map)
		readstats.printOtuMap(output.otumap_f);

	INFO("==== Done Post-processing (alignment statistics report) ====\n\n");
} // ~postProcess