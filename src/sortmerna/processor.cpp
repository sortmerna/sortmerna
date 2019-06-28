/**
 * FILE: processor.cpp
 * Created: Nov 26, 2017 Sun
 *
 * Callable objects designed to be run in threads
 *
 * @copyright 2016-19 Clarity Genomics BVBA
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
#include "read_control.hpp"
#include "writer.hpp"

// forward
void computeStats(Read & read, Readstats & readstats, Refstats & refstats, References & refs, Runopts & opts);

/* Runs in a thread. Pops reads from the Reads Queue */
void Processor::run()
{
	int countReads = 0;
	int countProcessed = 0;
	bool alreadyProcessed = false;
	std::stringstream ss;
	ss << "Processor " << id << " thread " << std::this_thread::get_id() << " started" << std::endl;
	std::cout << ss.str(); ss.str("");

	for (;;)
	{
		Read read = readQueue.pop(); // returns an empty read if queue is empty
		if (read.isEmpty && readQueue.getPushers() == 0)
		{
			break;
		}
		alreadyProcessed = (read.isRestored && read.lastIndex == index.index_num && read.lastPart == index.part);

		if (read.isEmpty || !read.isValid || alreadyProcessed) {
			if (alreadyProcessed) ++countProcessed;
			continue;
		}

		// search the forward and/or reverse strands depending on Run options
		int32_t num_strands = 0;
		//opts.forward = true; // TODO: this discards the possiblity of forward = false
		bool search_single_strand = opts.forward ^ opts.reverse; // search only a single strand
		if (search_single_strand)
			num_strands = 1; // only search the forward xor reverse strand
		else 
			num_strands = 2; // search both strands. The default when neither -F or -R were specified

		for (int32_t count = 0; count < num_strands; ++count)
		{
			if ((search_single_strand && opts.reverse) || count == 1)
			{
				if (!read.reversed)
					read.revIntStr();
			}
			// call 'paralleltraversal.cpp::alignmentCb'
			callback(opts, index, refs, output, readstats, refstats, read, search_single_strand || count == 1);
			//opts.forward = false;
			read.id_win_hits.clear(); // bug 46
		}

		if (read.isValid && !read.isEmpty) 
		{
			writeQueue.push(read);
		}

		countReads++;
	}
	writeQueue.decrPushers(); // signal this processor done adding
	writeQueue.notify(); // wake up writer waiting on queue.pop()

	ss << "Processor " << id << " thread " << std::this_thread::get_id() << " done. Processed " << countReads 
		<< " reads. Skipped already processed: " << countProcessed << " reads" << std::endl;
	std::cout << ss.str();
} // ~Processor::run

void PostProcessor::run()
{
	int countReads = 0;
	std::stringstream ss;

	ss << STAMP << "PostProcessor " << id << " thread " << std::this_thread::get_id() << " started" << std::endl;
	std::cout << ss.str();

	for (;;)
	{
		Read read = readQueue.pop(); // returns an empty read if queue is empty
		if (read.isEmpty)
		{ 
			if (readQueue.getPushers() == 0) 
				break; // queue is empty and no more pushers => end processing
			
			if (!read.isValid) 
				continue;
		}

		callback(read, readstats, refstats, refs, opts);
		++countReads;

		if (read.isValid && !read.isEmpty && !read.hit_denovo) 
		{
			writeQueue.push(read);
		}
	}
	writeQueue.decrPushers(); // signal this processor done adding
	writeQueue.notify(); // notify in case no Reads were ever pushed to the Write queue

	ss.str("");
	ss << STAMP << id << " thread " << std::this_thread::get_id() << " done. Processed [" << countReads << "] reads" << std::endl;
	std::cout << ss.str();
} // ~PostProcessor::run

void ReportProcessor::run()
{
	int countReads = 0;
	std::stringstream ss;

	ss << STAMP << "Report Processor " << id << " thread " << std::this_thread::get_id() << " started" << std::endl;
	std::cout << ss.str();

	int cap = opts.pairedin || opts.pairedout ? 2 : 1;
	std::vector<Read> reads;
	Read read;
	int i = 0;
	bool isDone = false;

	for (;!isDone;)
	{
		reads.clear();
		for (i = 0; i < cap; ++i)
		{
			read = readQueue.pop();  // returns an empty read if queue is empty
			reads.push_back(read);
			if (read.isEmpty)
			{
				if (readQueue.getPushers() == 0)
				{
					isDone = true;
					break;
				}
				if (!read.isValid)
					break;
			}
		}

		if (reads.back().isEmpty || !reads.back().isValid) continue;

		callback(reads, opts, refs, refstats, output);
		countReads+=i;
	}

	ss.str("");
	ss << STAMP << "Report Processor " << id << " thread " << std::this_thread::get_id() << " done. Processed [" << countReads << "] reads" << std::endl;
	std::cout << ss.str();
} // ~ReportProcessor::run

// called from main
void postProcess(Runopts & opts, Readstats & readstats, Output & output, KeyValueDatabase &kvdb)
{
	int N_READ_THREADS = opts.num_read_thread_pp;
	int N_PROC_THREADS = opts.num_proc_thread_pp; // opts.num_proc_threads
	int loopCount = 0; // counter of total number of processing iterations. TODO: no need here?
	std::stringstream ss;

	ss << STAMP << "\n\n==== Starting Post-processing routine (alignment statistics report) ====\n\n";
	std::cout << ss.str();

	ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS + opts.num_write_thread);
	ReadsQueue readQueue("read_queue", opts.queue_size_max, N_READ_THREADS); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", opts.queue_size_max, N_PROC_THREADS); // shared: Processor pushes, Writer pops
	bool indb = readstats.restoreFromDb(kvdb);

	if (indb) {
		ss.str("");
		ss << STAMP << "Restored Readstats from DB:\n    " << readstats.toString() << std::endl;
		std::cout << ss.str();
		ss.str("");
	}

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
				ss << STAMP << "Loading reference " << index_num 
					<< " part " << idx_part + 1 << "/" << refstats.num_index_parts[index_num] << "  ... ";
				std::cout << ss.str(); ss.str("");
				auto starts = std::chrono::high_resolution_clock::now(); // index loading start
				refs.load(index_num, idx_part, opts, refstats);
				std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
				ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << " sec]" << std::endl;
				std::cout << ss.str(); ss.str("");

				starts = std::chrono::high_resolution_clock::now(); // index processing starts

				for (int i = 0; i < N_READ_THREADS; ++i)
				{
					tpool.addJob(ReadControl(opts, readQueue, kvdb));
				}

				for (int i = 0; i < opts.num_write_thread; i++)
				{
					tpool.addJob(Writer("writer_" + std::to_string(i), writeQueue, kvdb, opts));
				}

				// add processor jobs
				for (int i = 0; i < N_PROC_THREADS; ++i)
				{
					tpool.addJob(PostProcessor("postproc_" + std::to_string(i), readQueue, writeQueue, opts, refs, readstats, refstats, computeStats));
				}
				++loopCount;
				tpool.waitAll(); // wait till processing is done on one index part
				refs.clear();
				readQueue.reset(N_READ_THREADS);
				writeQueue.reset(N_PROC_THREADS);

				elapsed = std::chrono::high_resolution_clock::now() - starts;
				ss << STAMP << "Done reference " << index_num << " Part: " << idx_part + 1
					<< " Time: " << std::setprecision(2) << std::fixed << elapsed.count() << " sec" << std::endl;
				std::cout << ss.str(); ss.str("");
			} // ~for(idx_part)
		} // ~for(index_num)

		ss.str("");
		ss << STAMP << "total_reads_denovo_clustering = " << readstats.total_reads_denovo_clustering << std::endl;
		std::cout << ss.str();

		readstats.set_is_total_reads_mapped_cov();
		readstats.is_stats_calc = true;
		readstats.store_to_db(kvdb); // store reads statistics computed by post-processor
	//} // ~if !readstats.stats_calc_done

	output.writeLog(opts, readstats);

	if (opts.otumapout)	readstats.printOtuMap(output.otumapFile);

	ss.str("");
	ss << STAMP << "\n\n==== Done Post-processing routine (alignment statistics report) ====\n\n";
	std::cout << ss.str();
} // ~postProcess