/* 
 * FILE: callbacks.cpp
 * Created: Jan 04, 2018 Thu
 *
 * @copyright 2016-19 Clarity Genomics BVBA
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


 // called on each read
void reportsJob(
	std::vector<Read> & reads, /* one or two (if paired) reads */
	Runopts & opts,
	References & refs,
	Refstats & refstats,
	Output & output
)
{
	// only needs one loop through all read, no reference file dependency
	if (opts.is_fastx && refs.num == 0 && refs.part == 0)
	{
		output.report_fasta(opts, reads);
	}

	// only needs one loop through all read, no reference file dependency
	if (opts.is_de_novo_otu && refs.num == 0 && refs.part == 0) {
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
void computeStats(Read & read, Readstats & readstats, Refstats & refstats, References & refs, Runopts & opts)
{
	// OTU-map: index of alignment holding maximum SW score
	uint32_t index_max_score = read.hits_align_info.max_index;
	if (read.is03) read.flip34();

	// loop all the alignments of this read
	for (uint32_t p = 0; p < read.hits_align_info.alignv.size(); ++p)
	{
		if (p == index_max_score)
		{
			// continue loop if the reference sequence in this alignment
			// belongs to the database section currently loaded into RAM
			if ((read.hits_align_info.alignv[p].index_num == refs.num) && (read.hits_align_info.alignv[p].part == refs.part))
			{
				// get the edit distance between reference and read
				uint32_t id = 0;
				uint32_t mismatches = 0;
				uint32_t gaps = 0;

				read.calcMismatchGapId(refs, p, mismatches, gaps, id);

				int32_t align_len = abs(read.hits_align_info.alignv[p].read_end1 + 1 - read.hits_align_info.alignv[p].read_begin1);
				int32_t total_pos = mismatches + gaps + id;
				std::stringstream ss;
				ss.precision(3);
				ss << (double)id / total_pos << ' ' << (double)align_len / read.hits_align_info.alignv[p].readlen;
				double align_id_round = 0.0;
				double align_cov_round = 0.0;
				ss >> align_id_round >> align_cov_round;

				// alignment with the highest SW score passed %id and %coverage thresholds
				if (align_id_round >= opts.align_id && align_cov_round >= opts.align_cov)
				{
					if (!readstats.is_total_reads_mapped_cov)
						++readstats.total_reads_mapped_cov; // if not already calculated.

					// TODO: this check is already performed during alignment (alignmentCb and compute_lis_alignment) 
					//       for (opts.num_alignments > -1)
					//  here for (opts.num_alignments == -1)
					// do not output read for de novo OTU construction (it passed the %id/coverage thresholds)
					if (opts.is_de_novo_otu) read.hit_denovo = false;

					// fill OTU map with highest-scoring alignment for the read
					if (opts.is_otu_map)
					{
						// reference sequence identifier for mapped read
						std::string refhead = refs.buffer[read.hits_align_info.alignv[p].ref_seq].header;
						std::string ref_seq_str = refhead.substr(0, refhead.find(' '));
						// left trim '>' or '@'
						ref_seq_str.erase(ref_seq_str.begin(),
							std::find_if(ref_seq_str.begin(), ref_seq_str.end(),
								[](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));

						// read identifier
						std::string read_seq_str = read.getSeqId();
						readstats.pushOtuMap(ref_seq_str, read_seq_str); // thread safe
					}
				} // ~if ID and Cov
			}//~if alignment at current database and index part loaded in RAM
			break; // no need to loop further after index_max_score was tested
		} // ~if p == index_max_score
	}//~for all alignments

	// only call once per read, on the last index/part
	if ( opts.is_de_novo_otu
		&& refs.num == opts.indexfiles.size() - 1 
		&& refs.part == refstats.num_index_parts[opts.indexfiles.size() - 1] -1 
		&& read.hit 
		&& read.hit_denovo )
	{
		++readstats.total_reads_denovo_clustering;
	}

} // ~computeStats