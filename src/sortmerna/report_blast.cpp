#include "report_blast.h"
#include "common.hpp"
#include "readfeed.hpp"
#include "options.hpp"
#include "read.hpp"
#include "references.hpp"
#include "refstats.hpp"

ReportBlast::ReportBlast(Runopts& opts)	: Report(opts) {}

ReportBlast::ReportBlast(Readfeed& readfeed, Runopts& opts)	: ReportBlast(opts)
{
	init(readfeed, opts);
}

void ReportBlast::init(Readfeed& readfeed, Runopts& opts) 
{
	fv.resize(readfeed.num_splits);
	fsv.resize(readfeed.num_splits);
	// WORKDIR/out/aligned_0_PID.blast
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + pid_str : "";
		fv[i] = opts.aligned_pfx.string() + sfx1 + sfx2 + ext;
		INFO("Testing file: ", std::filesystem::absolute(std::filesystem::path(fv[i])));
		fsv[i].open(fv[i]);
		fsv[i].close();
		if (!fsv[i]) {
			ERR("Failed stream on file ", fv[i]);
			exit(EXIT_FAILURE);
		}
	}
}

/**
 * called on each read => keep stream handle between calls
 */
void ReportBlast::append(int id, Read& read, References& refs, Refstats& refstats, Runopts& opts)
{
	const char MATCH = '|';
	const char MISMATCH = '*';
	const char INDEL = '-';

	uint32_t mid = 0; // matches id
	uint32_t mismatches = 0;
	uint32_t gaps = 0;
	char strandmark = '+';

	if (read.is03) read.flip34();

	// TODO: iterating all alignments for each reference part is an overhead. Alignments are pre-ordered, 
	//       so each new part corresponds to an index range of alignment vector. It's enough to loop 
	//       only that range.
	// iterate all alignments of the read
	for (int i = 0; i < read.alignment.alignv.size(); ++i)
	{
		if (read.alignment.alignv[i].index_num == refs.num
			&& read.alignment.alignv[i].part == refs.part)
		{
			// (λ*S - ln(K))/ln(2)
			uint32_t bitscore = (uint32_t)((float)((refstats.gumbel[refs.num].first)
				* (read.alignment.alignv[i].score1) - std::log(refstats.gumbel[refs.num].second)) / (float)std::log(2));

			// E = Kmn*exp(-λS)
			double evalue_score = (double)refstats.gumbel[refs.num].second
				* refstats.full_ref[refs.num]
				* refstats.full_read[refs.num]
				* std::exp(-refstats.gumbel[refs.num].first * read.alignment.alignv[i].score1);

			std::string refseq = refs.buffer[read.alignment.alignv[i].ref_num].sequence;
			std::string ref_id = refs.buffer[read.alignment.alignv[i].ref_num].id;

			if (read.alignment.alignv[i].strand)
				strandmark = '+';
			else
				strandmark = '-';

			if (read.alignment.alignv[i].strand == read.reversed) // XNOR
				read.revIntStr(); // reverse if necessary

			// Blast-like pairwise alignment (only for aligned reads)
			if (opts.blastFormat == BlastFormat::REGULAR)
			{
				fsv[id] << "Sequence ID: ";
				fsv[id] << ref_id; // print only start of the header till first space
				fsv[id] << std::endl;

				fsv[id] << "Query ID: ";
				fsv[id] << read.getSeqId();
				fsv[id] << std::endl;

				fsv[id] << "Score: " << read.alignment.alignv[i].score1 << " bits (" << bitscore << ")\t";
				fsv[id].precision(3);
				fsv[id] << "Expect: " << evalue_score << "\t";

				fsv[id] << "strand: " << strandmark << std::endl << std::endl;

				if (read.alignment.alignv[i].cigar.size() > 0)
				{
					uint32_t j, c = 0, left = 0, e = 0,
						qb = read.alignment.alignv[i].ref_begin1,
						pb = read.alignment.alignv[i].read_begin1;

					while (e < read.alignment.alignv[i].cigar.size() || left > 0)
					{
						int32_t count = 0;
						int32_t q = qb;
						int32_t p = pb;
						fsv[id] << "Target: ";
						fsv[id].width(8);
						fsv[id] << q + 1 << "    ";
						// process CIGAR
						for (c = e; c < read.alignment.alignv[i].cigar.size(); ++c)
						{
							// 4 Low bits encode a Letter: M | D | S
							uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
							// 28 High bits encode the number of occurencies e.g. 34
							uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 1) fsv[id] << INDEL; // mark indel
								else
								{
									fsv[id] << nt_map[(int)refseq[q]];
									++q;
								}
								++count;
								if (count == 60) goto step2;
							}
						}
					step2:
						fsv[id] << "    " << q << "\n";
						fsv[id].width(20);
						fsv[id] << " ";
						q = qb;
						count = 0;
						for (c = e; c < read.alignment.alignv[i].cigar.size(); ++c)
						{
							//uint32_t letter = 0xf & *(a->cigar + c);
							uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 0)
								{
									if ((char)nt_map[(int)refseq[q]] == (char)nt_map[(int)read.isequence[p]]) fsv[id] << MATCH; // mark match
									else fsv[id] << MISMATCH; // mark mismatch
									++q;
									++p;
								}
								else
								{
									fsv[id] << " ";
									if (letter == 1) ++p;
									else ++q;
								}
								++count;
								if (count == 60)
								{
									qb = q;
									goto step3;
								}
							}
						}
					step3:
						p = pb;
						fsv[id] << "\nQuery: ";
						fsv[id].width(9);
						fsv[id] << p + 1 << "    ";
						count = 0;
						for (c = e; c < read.alignment.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 2) fsv[id] << INDEL; // mark indel
								else
								{
									fsv[id] << nt_map[(int)read.isequence[p]];
									++p;
								}
								++count;
								if (count == 60)
								{
									pb = p;
									left = l - j - 1;
									e = (left == 0) ? (c + 1) : c;
									goto end;
								}
							}
						}
						e = c;
						left = 0;
					end:
						fsv[id] << "    " << p << "\n\n";
					}
				}
			}
			// Blast tabular m8 + optional columns for CIGAR and query coverage
			else if (opts.blastFormat == BlastFormat::TABULAR)
			{
				// (1) Query ID
				fsv[id] << read.getSeqId();

				// print null alignment for non-aligned read
				if (opts.is_print_all_reads && (read.alignment.alignv.size() == 0))
				{
					fsv[id] << "\t*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
					for (uint32_t l = 0; l < opts.blastops.size(); l++)
					{
						if (opts.blastops[l].compare("cigar") == 0)
							fsv[id] << "\t*";
						else if (opts.blastops[l].compare("qcov") == 0)
							fsv[id] << "\t0";
						else if (opts.blastops[l].compare("qstrand") == 0)
							fsv[id] << "\t*";
						fsv[id] << "\n";
					}
					return;
				}

				read.calcMismatchGapId(refs, i, mismatches, gaps, mid);
				int32_t total_pos = mismatches + gaps + mid;

				fsv[id] << "\t";
				// (2) Subject
				fsv[id] << ref_id << "\t";
				// (3) %id
				fsv[id].precision(3);
				fsv[id] << (double)mid / total_pos * 100 << "\t";
				// (4) alignment length
				fsv[id] << (read.alignment.alignv[i].read_end1 - read.alignment.alignv[i].read_begin1 + 1) << "\t";
				// (5) mismatches
				fsv[id] << mismatches << "\t";
				// (6) gap openings
				fsv[id] << gaps << "\t";
				// (7) q.start
				fsv[id] << read.alignment.alignv[i].read_begin1 + 1 << "\t";
				// (8) q.end
				fsv[id] << read.alignment.alignv[i].read_end1 + 1 << "\t";
				// (9) s.start
				fsv[id] << read.alignment.alignv[i].ref_begin1 + 1 << "\t";
				// (10) s.end
				fsv[id] << read.alignment.alignv[i].ref_end1 + 1 << "\t";
				// (11) e-value
				fsv[id] << evalue_score << "\t";
				// (12) bit score
				fsv[id] << bitscore;
				// OPTIONAL columns
				for (uint32_t l = 0; l < opts.blastops.size(); l++)
				{
					// output CIGAR string
					if (opts.blastops[l].compare("cigar") == 0)
					{
						fsv[id] << "\t";
						// masked region at beginning of alignment
						if (read.alignment.alignv[i].read_begin1 != 0) fsv[id] << read.alignment.alignv[i].read_begin1 << "S";
						for (int c = 0; c < read.alignment.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
							fsv[id] << length;
							if (letter == 0) fsv[id] << "M";
							else if (letter == 1) fsv[id] << "I";
							else fsv[id] << "D";
						}

						auto end_mask = read.sequence.length() - read.alignment.alignv[i].read_end1 - 1;
						// output the masked region at end of alignment
						if (end_mask > 0) fsv[id] << end_mask << "S";
					}
					// output % query coverage
					else if (opts.blastops[l].compare("qcov") == 0)
					{
						double coverage = (double)abs(read.alignment.alignv[i].read_end1 - read.alignment.alignv[i].read_begin1 + 1)
							/ read.alignment.alignv[i].readlen;

						fsv[id] << "\t";
						fsv[id].precision(3);
						fsv[id] << coverage * 100; // (double)align_len / readlen
					}
					// output strand
					else if (opts.blastops[l].compare("qstrand") == 0)
					{
						fsv[id] << "\t";
						fsv[id] << strandmark;
						//if (read.alignment.alignv[i].strand) blastout << "+";
						//else blastout << "-";
					}
				}
				fsv[id] << std::endl;
			}//~blast tabular m8
		}
	} // ~iterate all alignments
} // ~ ReportBlast::append