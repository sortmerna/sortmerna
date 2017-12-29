/**
* FILE: processor.cpp
* Created: Nov 26, 2017 Sun
*/

#include <iostream>
#include <sstream>

#include "processor.hpp"
#include "readsqueue.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "output.hpp"
#include "index.hpp"
#include "references.hpp"
#include "options.hpp"
#include "read.hpp"

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

		ss << "Processor: " << id << " Popped read id: " << read.id << " Index: " << read.lastIndex << " Part: " << read.lastPart << std::endl;
		std::cout << ss.str(); ss.str("");

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
			ss << "Processor: " << id << " Pushing read id: " << read.id << " Index: " << read.lastIndex << " Part: " << read.lastPart << std::endl;
			std::cout << ss.str(); ss.str("");
			writeQueue.push(read);
		}

		countReads++;
	}
	writeQueue.mDoneAdding();

	ss << "Processor " << id << " thread " << std::this_thread::get_id() << " done. Processed " << countReads
		<< " reads. Skipped (already processed) " << countProcessed << " reads\n";
	std::cout << ss.str(); ss.str("");
} // ~Processor::process

void ReportProcessor::run()
{
	int countReads = 0;
	std::stringstream ss;

	ss << "Report Processor " << id << " thread " << std::this_thread::get_id() << " started\n";
	std::cout << ss.str(); ss.str("");

	for (;!readQueue.isDone();)
	{
		Read read = readQueue.pop(); // returns an empty read if queue is empty

		if (read.isEmpty || !read.isValid)	continue;

		callback(opts, index, refs, output, readstats, refstats, read);

		countReads++;
	}

	ss << "Report Processor " << id << " thread " << std::this_thread::get_id() << " done. Processed " << countReads << " reads\n";
	std::cout << ss.str(); ss.str("");

} // ~ReportProcessor::process