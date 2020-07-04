/* 
 * FILE: callbacks.cpp
 * Created: Jan 04, 2018 Thu
 *
 * @copyright 2016-20 Clarity Genomics BVBA
 */

#include <string>
#include <cstdint>
#include <sstream>

#include "read.hpp"
#include "references.hpp"
#include "options.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "index.hpp"
#include "output.hpp"


 // called on each read or a pair of reads (if paired)
void reportsJob(std::vector<Read>& reads, Runopts& opts, References& refs, Refstats& refstats, Output& output)
{
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
} // ~reportsJob

/* 
 * Called for each index*index_part*read from PostProcessor::run
 *
 * Calculate:
 *     readstats.total_reads_mapped_cov 
 *     readstats.otu_map
 *     //read.hit_denovo see TODO in the function body
 */
void computeStats(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts)
{
	// OTU-map: index of alignment holding maximum SW score
	uint32_t index_max_score = read.alignment.max_index;
	if (read.is03) read.flip34();

	// populate OTU map
	if (opts.is_otu_map && read.is_id_cov) {
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

	// only call once per read, on the last index/part
	if ( opts.is_denovo_otu
		&& refs.num == opts.indexfiles.size() - 1 
		&& refs.part == refstats.num_index_parts[opts.indexfiles.size() - 1] -1 
		&& read.is_hit
		&& read.is_denovo)
	{
		++readstats.total_reads_denovo_clustering;
	}

} // ~computeStats