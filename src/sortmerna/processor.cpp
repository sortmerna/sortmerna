/**
 * FILE: processor.cpp
 * Created: Nov 26, 2017 Sun
 *
 * Callable objects designed to be run in threads
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
#include "reader.hpp"
#include "writer.hpp"

// forward
void computeStats(Read & read, Readstats & readstats, References & refs, Runopts & opts);
void writeLog(Runopts & opts, Readstats & readstats, Output & output);

void Processor::run()
{
	int countReads = 0;
	int countProcessed = 0;
	bool alreadyProcessed = false;
	std::stringstream ss;
	ss << "Processor " << id << " thread " << std::this_thread::get_id() << " started\n";
	std::cout << ss.str(); ss.str("");

	for (;!readQueue.isDone();)
	{
		Read read = readQueue.pop(); // returns an empty read if queue is empty
		alreadyProcessed = (read.isRestored && read.lastIndex == index.index_num && read.lastPart == index.part);

		//ss << "Processor: " << id << " Popped read id: " << read.id << " Index: " << read.lastIndex << " Part: " << read.lastPart << std::endl;
		//std::cout << ss.str(); ss.str("");

		if (read.isEmpty || !read.isValid || alreadyProcessed) {
			if (alreadyProcessed) ++countProcessed;
			continue;
		}

		// search the forward and/or reverse strands depending on Run options
		int32_t strandCount = 0;
		opts.forward = true; // TODO: this discards the possiblity of forward = false
		if (opts.forward ^ opts.reverse)
			strandCount = 1; // only search the forward xor reverse strand
		else 
			strandCount = 2; // search both strands. The default when neither -F or -R were specified

		for (int32_t strand = 0; strand < strandCount; strand++)
		{
			if (!opts.forward && !read.reversed)
				read.revIntStr(); // reverse-complement the sequence
			callback(opts, index, refs, output, readstats, refstats, read);
			opts.forward = false;
		}

		if (read.isValid && !read.isEmpty) {
			//ss << "Processor: " << id << " Pushing read id: " << read.id << " Index: " << read.lastIndex << " Part: " << read.lastPart << std::endl;
			//std::cout << ss.str(); ss.str("");
			writeQueue.push(read);
		}

		countReads++;
	}
	writeQueue.mDoneAdding();

	ss << "Processor " << id << " thread " << std::this_thread::get_id() << " done. Processed " << countReads
		<< " reads. Skipped (already processed) " << countProcessed << " reads\n";
	std::cout << ss.str(); ss.str("");
} // ~Processor::run

void PostProcessor::run()
{
	int countReads = 0;
	std::stringstream ss;

	ss << "PostProcessor " << id << " thread " << std::this_thread::get_id() << " started\n";
	std::cout << ss.str(); ss.str("");

	for (;!readQueue.isDone();)
	{
		Read read = readQueue.pop(); // returns an empty read if queue is empty

		if (read.isEmpty || !read.isValid)	continue;

		callback(read, readstats, refs, opts);
		++countReads;

		if (read.isValid && !read.isEmpty && !read.hit_denovo) 
		{
			writeQueue.push(read);
		}
	}
	writeQueue.mDoneAdding();

	ss << "PostProcessor " << id << " thread " << std::this_thread::get_id() << " done. Processed " << countReads << " reads\n";
	std::cout << ss.str(); ss.str("");
} // ~PostProcessor::run

void ReportProcessor::run()
{
	int countReads = 0;
	std::stringstream ss;

	ss << "Report Processor " << id << " thread " << std::this_thread::get_id() << " started\n";
	std::cout << ss.str(); ss.str("");
	int cap = opts.pairedin || opts.pairedout ? 2 : 1;
	std::vector<Read> reads;
	int i = 0;

	for (;!readQueue.isDone();)
	{
		for (i = 0; i < cap; ++i)
		{
			reads.push_back(readQueue.pop()); // returns an empty read if queue is empty
			if (reads[i].isEmpty || !reads[i].isValid) break;
		}

		if (reads.back().isEmpty || !reads.back().isValid)	continue;

		callback(reads, opts, refs, refstats, output);
		reads.clear();
		countReads+=i;
	}

	ss << "Report Processor " << id << " thread " << std::this_thread::get_id() << " done. Processed " << countReads << " reads\n";
	std::cout << ss.str(); ss.str("");

} // ~ReportProcessor::run

// called from main
void postProcess(Runopts & opts, Readstats & readstats, Output & output)
{
	int N_READ_THREADS = opts.num_read_thread_pp;
	int N_PROC_THREADS = opts.num_proc_thread_pp; // opts.num_proc_threads
	int loopCount = 0; // counter of total number of processing iterations. TODO: no need here?
	std::stringstream ss;

	ss << "\tpostProcess Thread: " << std::this_thread::get_id() << std::endl;
	std::cout << ss.str(); ss.str("");

	ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS);
	KeyValueDatabase kvdb(opts.kvdbPath);
	ReadsQueue readQueue("read_queue", QUEUE_SIZE_MAX, N_READ_THREADS); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", QUEUE_SIZE_MAX, N_PROC_THREADS); // shared: Processor pushes, Writer pops
	readstats.restoreFromDb(kvdb);
	Refstats refstats(opts, readstats);
	References refs;

	// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); ++index_num)
	{
		// iterate parts of reference files
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
		{
			ss << "\tpostProcess: Loading reference " << index_num << " part " << idx_part + 1 << "/" << refstats.num_index_parts[index_num] << "  ... ";
			std::cout << ss.str(); ss.str("");
			auto starts = std::chrono::high_resolution_clock::now(); // index loading start
			refs.load(index_num, idx_part, opts, refstats);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
			ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << " sec]\n";
			std::cout << ss.str(); ss.str("");

			starts = std::chrono::high_resolution_clock::now(); // index processing starts

			for (int i = 0; i < N_READ_THREADS; ++i)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), opts, readQueue, kvdb, loopCount));
			}

			for (int i = 0; i < opts.num_write_thread; i++)
			{
				tpool.addJob(Writer("writer_" + std::to_string(i), writeQueue, kvdb));
			}

			// add processor jobs
			for (int i = 0; i < N_PROC_THREADS; ++i)
			{
				tpool.addJob(PostProcessor("postproc_" + std::to_string(i), readQueue, writeQueue, opts, refs, readstats, computeStats));
			}
			++loopCount;
			tpool.waitAll(); // wait till processing is done on one index part
			refs.clear();
			readQueue.reset(N_READ_THREADS);
			writeQueue.reset(N_PROC_THREADS);

			elapsed = std::chrono::high_resolution_clock::now() - starts;
			ss << "    Done reference " << index_num << " Part: " << idx_part + 1
				<< " Time: " << std::setprecision(2) << std::fixed << elapsed.count() << " sec\n";
			std::cout << ss.str(); ss.str("");
		} // ~for(idx_part)
	} // ~for(index_num)

	ss << "readstats.total_reads_denovo_clustering: " << readstats.total_reads_denovo_clustering << std::endl;
	std::cout << ss.str(); ss.str("");

	writeLog(opts, readstats, output);

	if (opts.otumapout)	readstats.printOtuMap(output.otumapFile);

	kvdb.put("Readstats", readstats.toString()); // store statistics computed by post-processor

	std::cout << "\tpostProcess: Done \n";
} // ~postProcess

void writeLog(Runopts & opts, Readstats & readstats, Output & output)
{
	output.openfiles(opts);

	if (opts.samout) output.writeSamHeader(opts);

	// output total number of reads
	output.logstream << " Results:\n";
	output.logstream << "    Total reads = " << readstats.number_total_read << "\n";
	if (opts.de_novo_otu)
	{
		// all reads that have read::hit_denovo == true
		output.logstream << "    Total reads for de novo clustering = " << readstats.total_reads_denovo_clustering << "\n";
	}
	// output total non-rrna + rrna reads
	output.logstream << std::setprecision(2) << std::fixed;
	output.logstream << "    Total reads passing E-value threshold = " << readstats.total_reads_mapped.load()
		<< " (" << (float)((float)readstats.total_reads_mapped.load() / (float)readstats.number_total_read) * 100 << ")\n";
	output.logstream << "    Total reads failing E-value threshold = "
		<< readstats.number_total_read - readstats.total_reads_mapped.load()
		<< " (" << (1 - ((float)((float)readstats.total_reads_mapped.load() / (float)readstats.number_total_read))) * 100 << ")\n";
	output.logstream << "    Minimum read length = " << readstats.min_read_len << "\n";
	output.logstream << "    Maximum read length = " << readstats.max_read_len << "\n";
	output.logstream << "    Mean read length    = " << readstats.full_read_main / readstats.number_total_read << "\n";

	output.logstream << " By database:\n";

	// output stats by database
	for (uint32_t index_num = 0; index_num < opts.indexfiles.size(); index_num++)
	{
		output.logstream << "    " << opts.indexfiles[index_num].first << "\t\t"
			<< (float)((float)readstats.reads_matched_per_db[index_num] / (float)readstats.number_total_read) * 100 << "\n";
	}

	if (opts.otumapout)
	{
		output.logstream << " Total reads passing %%id and %%coverage thresholds = " << readstats.total_reads_mapped_cov.load() << "\n";
		output.logstream << " Total OTUs = " << readstats.otu_map.size() << "\n";
	}
	time_t q = time(0);
	struct tm * now = localtime(&q);
	output.logstream << "\n " << asctime(now) << "\n";
} // ~writeLog