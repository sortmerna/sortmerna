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
#include <functional> // std::ref

#include "processor.hpp"
#include "readsqueue.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "output.hpp"
#include "index.hpp"
#include "references.hpp"
#include "options.hpp"
#include "read.hpp"
//#include "ThreadPool.hpp"
#include "readfeed.hpp"
#include "summary.hpp"

// forward
void writeSummary3(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts);
void traverse(Runopts& opts, Index& index, References& refs, Readstats& readstats, Refstats& refstats, Read& read, bool isLastStrand);

/*
 * Called for each index*index_part*read
 *
 * populate 'readstats.otu_map'
 * count 'readstats.total_reads_denovo_clustering'
 */
void writeSummary3(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts)
{
	//uint32_t index_max_score = read.alignment.max_index; // index of alignment holding maximum SW score

	if (opts.is_otu_map || opts.is_denovo_otu) {
		if (read.is03) read.flip34();
	}

	// populate OTU map
	if (opts.is_otu_map) {
		if (read.is_id && read.is_cov) {
			// reference sequence identifier for mapped read
			std::string refhead = refs.buffer[read.alignment.alignv[read.alignment.max_index].ref_num].header;
			std::string ref_seq_str = refhead.substr(0, refhead.find(' '));
			// left trim '>' or '@'
			ref_seq_str.erase(ref_seq_str.begin(),
				std::find_if(ref_seq_str.begin(), ref_seq_str.end(),
					[](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));

			// read identifier
			std::string read_seq_str = read.getSeqId();
			readstats.pushOtuMap(ref_seq_str, read_seq_str); // thread safe
		}
	}

	// only call once per read, on the last index/part
	if (opts.is_denovo_otu) {
		if (refs.num == opts.indexfiles.size() - 1
			&& refs.part == refstats.num_index_parts[opts.indexfiles.size() - 1] - 1
			&& read.is_hit
			&& read.is_denovo)
		{
			++readstats.total_reads_denovo_clustering;
		}
	}
} // ~postProcess3

/*
  runs in a thread
*/
void writeSummary2(int id, Readfeed& readfeed, Runopts& opts, References& refs, Readstats& readstats, Refstats& refstats, KeyValueDatabase& kvdb)
{
	unsigned countReads = 0;
	unsigned count_reads_aligned = 0;
	std::string readstr;

	INFO("PostProcessor: ", id, " thread: ", std::this_thread::get_id(), " started");

	for (;readfeed.next(id, readstr);)
	{
		{
			Read read(readstr);
			read.init(opts);
			read.load_db(kvdb);

			if (!read.isValid)
				continue;

			writeSummary3(read, readstats, refstats, refs, opts);
			readstr.resize(0);
			++countReads;
			if (read.is_hit) ++count_reads_aligned;
		}
	}

	INFO("PostProcessor: ", id, " thread: ", std::this_thread::get_id(),
		" done. Processed reads: ", countReads, ". count_reads_aligned: ", count_reads_aligned);
}

// called from main
void writeSummary(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	INFO("==== Starting summary (alignment statistics report) ====\n\n");
	Summary summary;
	bool indb = readstats.restoreFromDb(kvdb);
	if (indb) 
		INFO("Restored Readstats from DB:\n    ", readstats.toString());

	readstats.total_reads_denovo_clustering = 0; // TODO: to prevent incrementing the stored value. Change this if ever using 'stats_calc_done"

	//if (!readstats.stats_calc_done)
	//{
		Refstats refstats(opts, readstats);

		// this part is only necessary for OTU map and/or deNovo clustering
		if (opts.is_otu_map || opts.is_denovo_otu) {
			//ReadsQueue read_queue("queue_1", opts.queue_size_max, readstats.all_reads_count);
			int numThreads = 0;
			if (opts.feed_type == FEED_TYPE::LOCKLESS)
			{
				numThreads = opts.num_read_thread_pp + opts.num_proc_thread_pp;
				INFO("using total threads: ", numThreads, " including Read threads: ", opts.num_read_thread_pp, " Processor threads: ", opts.num_proc_thread_pp);
			}
			else {
				numThreads = opts.num_proc_thread_pp;
				INFO("Using total threads: ", numThreads);
			}

			std::vector<std::thread> tpool;
			tpool.reserve(numThreads);

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
					tpool.emplace_back(std::thread(writeSummary2, 0, std::ref(readfeed), std::ref(opts), std::ref(refs),
						std::ref(readstats), std::ref(refstats), std::ref(kvdb)));

					// wait till processing is done on one index part
					//tpool.waitAll(); 
					for (auto i = 0; i < tpool.size(); ++i) {
						tpool[i].join();
					}

					refs.unload();
					//read_queue.reset();

					elapsed = std::chrono::high_resolution_clock::now() - starts;
					INFO_MEM("Done reference ", idx, " Part: ", ipart + 1, " Elapsed sec: ", elapsed.count());
				} // ~for(ipart)
			} // ~for(idx)

			INFO("total_reads_denovo_clustering = ", readstats.total_reads_denovo_clustering);
		} // ~if opts.is_otu_map || opts.is_denovo_otu

		readstats.set_is_total_mapped_sw_id_cov();
		readstats.is_stats_calc = true;
		readstats.store_to_db(kvdb); // store reads statistics computed by post-processor
	//} // ~if !readstats.stats_calc_done

	summary.write(opts, refstats, readstats);

	if (opts.is_otu_map) {
		// OTU map output file  WORKDIR/out/aligned_otus.txt
		std::string f_otumap;
		std::ofstream otumap;
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + std::to_string(getpid());
		}
		sfx += "_otus.txt";
		f_otumap = opts.aligned_pfx.string() + sfx;
		readstats.printOtuMap(f_otumap);
	}

	INFO("==== Done Post-processing (alignment statistics report) ====\n\n");
} // ~writeSummary

// 20201004 moved here from callbacks.cpp
void report(int id, Readfeed& readfeed, 
	            Runopts& opts, 
	            References& refs, 
	            Refstats& refstats, 
	            Output& output, 
	            KeyValueDatabase& kvdb)
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
			if (readfeed.next(id, readstr))
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
				//job(reads, opts, refs, refstats, output);

	// only needs one loop through all read, no reference file dependency
				if (opts.is_fastx && refs.num == 0 && refs.part == 0)
				{
					output.report_fasta(opts, reads);
				}

				// only needs one loop through all read, no reference file dependency
				if (opts.is_denovo_otu && refs.num == 0 && refs.part == 0) {
					output.report_denovo(opts, reads);
				}

				for (Read read : reads)
				{
					if (opts.is_blast)
					{
						output.report_blast(opts, refstats, refs, read);
					}

					if (opts.is_sam)
					{
						output.report_sam(opts, refs, read);
					}
				} // ~for reads
			}
		}
	} // ~for

	INFO_MEM("Report Processor: ", id, " thread: ", std::this_thread::get_id(), " done. Processed reads: ", countReads, " Invalid reads: ", num_invalid);
} // ~report


// called from main. generateReports -> reportsJob
void writeReports(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	INFO("=== Report generation starts. Thread: ", std::this_thread::get_id(), " ===\n");

	int nthreads = 0;
	if (readfeed.type == FEED_TYPE::SPLIT_READS) {
		nthreads = opts.num_proc_thread;
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

	output.openfiles(opts);
	if (opts.is_sam) output.writeSamHeader(opts);

	// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t ref_idx = 0; ref_idx < (uint16_t)opts.indexfiles.size(); ++ref_idx)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[ref_idx]; ++idx_part)
		{
			INFO("Loading reference ", ref_idx, " part ", idx_part + 1, "/", refstats.num_index_parts[ref_idx], "  ... ");

			auto starts = std::chrono::high_resolution_clock::now();

			refs.load(ref_idx, idx_part, opts, refstats);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts; // ~20 sec Debug/Win
			INFO("done. Elapsed sec [", elapsed.count());

			starts = std::chrono::high_resolution_clock::now(); // index processing starts

			// start processor
			if (opts.feed_type == FEED_TYPE::SPLIT_READS) {
				for (int i = 0; i < tpool.size(); ++i) {
					tpool.emplace_back(std::thread(report, i, std::ref(readfeed), std::ref(opts), 
						std::ref(refs), std::ref(refstats), std::ref(output), std::ref(kvdb)));
				}
			}
			// wait till processing is done
			for (auto i = 0; i < tpool.size(); ++i) {
				tpool[i].join();
			}

			refs.unload();
			//read_queue.reset();

			elapsed = std::chrono::high_resolution_clock::now() - starts; // index processing done
			INFO("Done reference ", ref_idx, " Part: ", idx_part + 1, " Elapsed sec: ", elapsed.count());

			if (!opts.is_blast && !opts.is_sam)	break;;
		} // ~for(idx_part)
	} // ~for(ref_idx)

	INFO("=== Done Reports generation ===\n");
} // ~generateReports

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