/**
 * FILE: processor.cpp
 * Created: Nov 26, 2017 Sun
 *
 * performs the alignment
 *
 * @copyright 2016-20 Clarity Genomics BVBA
 */

#include <chrono>

#include "processor.hpp"
#include "read.hpp"
#include "readfeed.hpp"
#include "index.hpp"
#include "references.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "options.hpp"
//#include "readsqueue.hpp"

// forward
void traverse(Runopts& opts, Index& index, References& refs, Readstats& readstats, Refstats& refstats, Read& read, bool isLastStrand);

/*
  runs in a thread.  align -> align2
  @param  id
*/
void align2(int id, Readfeed& readfeed, Readstats& readstats, Index& index, References& refs, Refstats& refstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	unsigned num_all = 0; // all reads this processor sees
	unsigned num_skipped = 0; // reads already processed i.e. results found in Database
	unsigned num_hit = 0; // count of reads with read.hit = true found by a single thread - just for logging
	std::string readstr;

	INFO("Processor ", id, " thread ", std::this_thread::get_id(), " started");
	auto idx = id * readfeed.num_orig_files;
	for (auto incr = 1; readfeed.next(idx, readstr);)
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
				
				traverse(opts, index, refs, readstats, refstats, read, search_single_strand || count == 1); // 'paralleltraversal.cpp'
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

		// switch FWD-REV if two files are processed
		if (readfeed.is_two_files) {
			idx += incr;
			incr = incr > 0 ? -1 : 1;
		}
	} // ~while there are reads

	INFO("Processor ", id, " thread ", std::this_thread::get_id(), " done. Processed ", num_all,
		" reads. Skipped already processed: ", num_skipped, " reads", " Aligned reads (passing E-value): ", num_hit);
} // ~align2

// called from main
void align(Readfeed& readfeed, Readstats& readstats, Index& index, KeyValueDatabase& kvdb, Runopts& opts)
{
	INFO("==== Starting alignment ====");

	unsigned int numCores = std::thread::hardware_concurrency(); // find number of CPU cores
	INFO("Number of cores: ", numCores);

	// Init thread pool with the given number of threads
	int numProcThread = 0;
	numProcThread = opts.num_proc_thread; // '-thread'

	// calculate the number of threads to use
	int numThreads = 0;
	if (opts.feed_type == FEED_TYPE::LOCKLESS)
	{
		numThreads = opts.num_read_thread + numProcThread;
		INFO("using total threads: ", numThreads, " including Read threads: ", opts.num_read_thread, " Processor threads: ", numProcThread);
		//ThreadPool tpool(numThreads);
		//ReadsQueue read_queue("queue_1", opts.queue_size_max, readstats.all_reads_count, numProcThread);
	}
	else {
		numThreads = numProcThread;
		INFO("Using number of Processor threads: ", numProcThread);
		readfeed.init_reading(); // prepare readfeed
	}
	std::vector<std::thread> tpool;
	tpool.reserve(numThreads);

	Refstats refstats(opts, readstats);
	References refs;

	int loopCount = 0; // counter of total number of processing iterations

	// perform alignment
	auto starts = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;

	// loop through every index passed to option '--ref'
	for (size_t index_num = 0; index_num < opts.indexfiles.size(); ++index_num)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
		{
			// load index
			INFO("Loading index: ", index_num, " part: ", idx_part + 1, "/", refstats.num_index_parts[index_num], " Memory KB: ", (get_memory() >> 10), " ... ");
			starts = std::chrono::high_resolution_clock::now();
			index.load(index_num, idx_part, opts.indexfiles, refstats);
			readstats.short_reads_num.store(0, std::memory_order_relaxed); // reset the short reads counter
			elapsed = std::chrono::high_resolution_clock::now() - starts; // ~20 sec Debug/Win
			INFO_MEM("done [", elapsed.count(), "] sec");

			// load references
			INFO("Loading references ...");
			starts = std::chrono::high_resolution_clock::now();
			refs.load(index_num, idx_part, opts, refstats);
			elapsed = std::chrono::high_resolution_clock::now() - starts; // ~20 sec Debug/Win
			INFO_MEM("done [", elapsed.count(), "] sec.");

			starts = std::chrono::high_resolution_clock::now();

			// add Readfeed job if necessary
			if (opts.feed_type == FEED_TYPE::LOCKLESS)
			{
				//tpool.addJob(f_readfeed_run);
			}

			// add Processor jobs
			for (int i = 0; i < numProcThread; i++)
			{
				tpool.emplace_back(std::thread(align2, i, std::ref(readfeed), std::ref(readstats), std::ref(index),
					std::ref(refs), std::ref(refstats), std::ref(kvdb), std::ref(opts)));
			}
			for (auto i = 0; i < tpool.size(); ++i) {
				tpool[i].join();
			}

			++loopCount;

			elapsed = std::chrono::high_resolution_clock::now() - starts;
			//INFO_MEM("Done index ", index_num, " Part: ", idx_part + 1, " Queue size: ", read_queue.queue.size_approx(), " Time: ", elapsed.count())

			index.unload();
			refs.unload();
			INFO_MEM("Index and References unloaded.");
			tpool.clear();
			// rewind for the next index
			readfeed.rewind_in();
			readfeed.init_vzlib_in();
			//read_queue.reset();
		} // ~for(idx_part)
	} // ~for(index_num)

	INFO("==== Done alignment ====\n");

	// store readstats calculated in alignment
	readstats.set_is_total_mapped_sw_id_cov();
	readstats.store_to_db(kvdb);
} // ~align