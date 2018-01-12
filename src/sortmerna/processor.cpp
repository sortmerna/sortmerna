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

// forward
void computeStats(Read & read, Readstats & readstats, References & refs, Runopts & opts);
void writeLog(Runopts & opts, Readstats & readstats);

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
				read.revIntStr(); // reverse the sequence
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

	for (;!readQueue.isDone(); countReads++)
	{
		Read read = readQueue.pop(); // returns an empty read if queue is empty

		if (read.isEmpty || !read.isValid)	continue;

		callback(read, readstats, refs, opts);
	}

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

void runPostProcessor(Runopts & opts)
{
	int N_READ_THREADS = 1;
	int N_PROC_THREADS = 1;
	int loopCount = 0; // counter of total number of processing iterations. TODO: no need here?
	std::stringstream ss;

	ss << "runPostProcessor Thread: " << std::this_thread::get_id() << std::endl;
	std::cout << ss.str(); ss.str("");

	ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS);
	KeyValueDatabase kvdb(opts.kvdbPath);
	ReadsQueue readQueue("read_queue", QUEUE_SIZE_MAX, N_READ_THREADS); // shared: Processor pops, Reader pushes
	Readstats readstats(opts);
	readstats.restoreFromDb(kvdb);
	Refstats refstats(opts, readstats);
	References refs;

	// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); ++index_num)
	{
		// iterate parts of reference files
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
		{
			ss << "runPostProcessor: Loading reference part " << idx_part + 1 << "/" << refstats.num_index_parts[index_num] << "  ... ";
			std::cout << ss.str(); ss.str("");
			auto t = std::chrono::high_resolution_clock::now();
			refs.load(index_num, idx_part, opts, refstats);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t;
			ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << " sec]\n";
			std::cout << ss.str(); ss.str("");

			for (int i = 0; i < N_READ_THREADS; ++i)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), opts, readQueue, kvdb, loopCount));
			}

			// add processor jobs
			for (int i = 0; i < N_PROC_THREADS; ++i)
			{
				tpool.addJob(PostProcessor("postproc_" + std::to_string(i), readQueue, opts, refs, readstats, computeStats));
			}
			++loopCount;
			tpool.waitAll(); // wait till processing is done on one index part
			refs.clear();
			readQueue.reset(N_READ_THREADS);
		} // ~for(idx_part)
	} // ~for(index_num)
	writeLog(opts, readstats);
	kvdb.put("Readstats", readstats.toString()); // store statistics computed by post-processor
	std::cout << "runPostProcessor: Done \n";
} // ~runPostProcessor

void writeLog(Runopts & opts, Readstats & readstats)
{
	Output output(opts, readstats);
	output.openfiles(opts);
	if (opts.samout) output.writeSamHeader(opts);
	output.logstream.open(output.logfile, std::ofstream::binary | std::ofstream::app);

	// output total number of reads
	output.logstream << " Results:\n";
	output.logstream << "    Total reads = " << readstats.number_total_read << "\n";
	if (opts.de_novo_otu)
	{
		// total_reads_denovo_clustering = sum of all reads that have read::hit_denovo == true
		// either query DB or store in Readstats::total_reads_denovo_clustering
		output.logstream << "    Total reads for de novo clustering = " << readstats.total_reads_denovo_clustering << "\n";
	}
	// output total non-rrna + rrna reads
	output.logstream << std::setprecision(2) << std::fixed;
	output.logstream << "    Total reads passing E-value threshold = " << readstats.total_reads_mapped
		<< " (" << (float)((float)readstats.total_reads_mapped / (float)readstats.number_total_read) * 100 << ")\n";
	output.logstream << "    Total reads failing E-value threshold = "
		<< readstats.number_total_read - readstats.total_reads_mapped
		<< " (" << (1 - ((float)((float)readstats.total_reads_mapped / (float)readstats.number_total_read))) * 100 << ")\n";
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
		output.logstream << " Total reads passing %%id and %%coverage thresholds = " << readstats.total_reads_mapped_cov << "\n";
		output.logstream << " Total OTUs = " << readstats.otu_map.size() << "\n";
	}
	time_t q = time(0);
	struct tm * now = localtime(&q);
	output.logstream << "\n " << asctime(now) << "\n";
	output.logstream.close();
} // ~writeLog