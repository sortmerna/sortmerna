/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/
#include <cmath> // std::log, std::exp

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
	openfw();
}

void ReportBlast::init(Readfeed& readfeed, Runopts& opts) 
{
	fv.resize(readfeed.num_splits);
	fsv.resize(readfeed.num_splits);
	is_zip = readfeed.orig_files[0].isZip;
	// WORKDIR/out/aligned_0_PID.blast
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + pid_str : "";
		std::string gz = is_zip ? ".gz" : "";
		fv[i] = opts.aligned_pfx.string() + sfx1 + sfx2 + ext + gz;
		openfw(i);
	}
	if (is_zip) init_zip();
}

/**
 * called on each read => keep stream handle between calls
 */
void ReportBlast::append(int id, Read& read, References& refs, Refstats& refstats, Runopts& opts)
{
	const char MATCH = '|';
	const char MISMATCH = '*';
	const char INDEL = '-';
	char strandmark = '+';
	std::stringstream ss;

	if (read.is03) read.flip34();

	// TODO: iterating all alignments for each reference part is an overhead. Alignments are pre-ordered, 
	//       so each new part corresponds to an index range of alignment vector. It's enough to loop 
	//       only that range.
	// iterate all alignments of the read
	for (auto const& align: read.alignment.alignv)
	{
		if (align.index_num == refs.num
			&& align.part == refs.part)
		{
			// (λ*S - ln(K))/ln(2)
			uint32_t bitscore = (uint32_t)((float)((refstats.gumbel[refs.num].first)
				* (align.score1) - std::log(refstats.gumbel[refs.num].second)) / (float)std::log(2));

			// E = Kmn*exp(-λS)
			double evalue_score = (double)refstats.gumbel[refs.num].second
				* refstats.full_ref[refs.num]
				* refstats.full_read[refs.num]
				* std::exp(-refstats.gumbel[refs.num].first * align.score1);

			std::string refseq = refs.buffer[align.ref_num].sequence;
			std::string ref_id = refs.buffer[align.ref_num].id;

			strandmark = align.strand ? '+' : '-';

			if (align.strand == read.reversed) // XNOR
				read.revIntStr(); // reverse if necessary

			// Blast-like pairwise alignment (only for aligned reads)
			if (opts.blastFormat == BlastFormat::REGULAR)
			{
				ss << "Sequence ID: " << ref_id << std::endl; // print only start of the header till first space
				ss << "Query ID: " << read.getSeqId() << std::endl;

				ss << "Score: " << align.score1 << " bits (" << bitscore << ")\t";
				ss.precision(3);
				ss << "Expect: " << evalue_score << "\t";

				ss << "strand: " << strandmark << std::endl << std::endl;

				if (align.cigar.size() > 0)
				{
					uint32_t j, c = 0, left = 0, e = 0,
						qb = align.ref_begin1,
						pb = align.read_begin1;

					while (e < align.cigar.size() || left > 0)
					{
						int32_t count = 0;
						int32_t q = qb;
						int32_t p = pb;
						ss << "Target: ";
						ss.width(8);
						ss << q + 1 << "    ";
						// process CIGAR
						for (c = e; c < align.cigar.size(); ++c)
						{
							// 4 Low bits encode a Letter: M | D | S
							uint32_t letter = 0xf & align.cigar[c];
							// 28 High bits encode the number of occurencies e.g. 34
							uint32_t length = (0xfffffff0 & align.cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 1) ss << INDEL; // mark indel
								else
								{
									ss << nt_map[(int)refseq[q]];
									++q;
								}
								++count;
								if (count == 60) goto step2;
							}
						}
					step2:
						ss << "    " << q << "\n";
						ss.width(20);
						ss << " ";
						q = qb;
						count = 0;
						for (c = e; c < align.cigar.size(); ++c)
						{
							//uint32_t letter = 0xf & *(a->cigar + c);
							uint32_t letter = 0xf & align.cigar[c];
							uint32_t length = (0xfffffff0 & align.cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 0)
								{
									if ((char)nt_map[(int)refseq[q]] == (char)nt_map[(int)read.isequence[p]]) ss << MATCH; // mark match
									else ss << MISMATCH; // mark mismatch
									++q;
									++p;
								}
								else
								{
									ss << " ";
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
						ss << "\nQuery: ";
						ss.width(9);
						ss << p + 1 << "    ";
						count = 0;
						for (c = e; c < align.cigar.size(); ++c)
						{
							uint32_t letter = 0xf & align.cigar[c];
							uint32_t length = (0xfffffff0 & align.cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 2) ss << INDEL; // mark indel
								else
								{
									ss << nt_map[(int)read.isequence[p]];
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
						ss << "    " << p << "\n\n";
					}
				}
			}
			// Blast tabular m8 + optional columns for CIGAR and query coverage
			else if (opts.blastFormat == BlastFormat::TABULAR)
			{
				// (1) Query ID
				ss << read.getSeqId();

				// print null alignment for non-aligned read
				if (opts.is_print_all_reads && (read.alignment.alignv.size() == 0))
				{
					ss << "\t*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
					for (uint32_t l = 0; l < opts.blastops.size(); l++)
					{
						if (opts.blastops[l].compare("cigar") == 0)
							ss << "\t*";
						else if (opts.blastops[l].compare("qcov") == 0)
							ss << "\t0";
						else if (opts.blastops[l].compare("qstrand") == 0)
							ss << "\t*";
						ss << "\n";
					}
					return;
				}

				auto miss_gap_match = read.calc_miss_gap_match(refs, align);

				ss << "\t";
				// (2) Subject
				ss << ref_id << "\t";
				// (3) %id
				ss.precision(3);
				ss << std::get<3>(miss_gap_match) * 100 << "\t";
				// (4) alignment length
				ss << (align.read_end1 - align.read_begin1 + 1) << "\t";
				// (5) mismatches
				ss << std::get<0>(miss_gap_match) << "\t";
				// (6) gap openings
				ss << std::get<1>(miss_gap_match) << "\t";
				// (7) q.start
				ss << align.read_begin1 + 1 << "\t";
				// (8) q.end
				ss << align.read_end1 + 1 << "\t";
				// (9) s.start
				ss << align.ref_begin1 + 1 << "\t";
				// (10) s.end
				ss << align.ref_end1 + 1 << "\t";
				// (11) e-value
				ss << evalue_score << "\t";
				// (12) bit score
				ss << bitscore;
				// OPTIONAL columns: CIGAR, %COV, strand
				for (uint32_t l = 0; l < opts.blastops.size(); l++)
				{
					if (opts.blastops[l].compare("cigar") == 0)
					{
						// output CIGAR string
						ss << "\t";
						// masked region at beginning of alignment
						if (align.read_begin1 != 0) ss << align.read_begin1 << "S";
						for (int c = 0; c < align.cigar.size(); ++c)
						{
							uint32_t letter = 0xf & align.cigar[c];
							uint32_t length = (0xfffffff0 & align.cigar[c]) >> 4;
							ss << length;
							if (letter == 0) ss << "M";
							else if (letter == 1) ss << "I";
							else ss << "D";
						}

						auto end_mask = read.sequence.length() - align.read_end1 - 1;
						// output the masked region at end of alignment
						if (end_mask > 0) ss << end_mask << "S";
					}
					else if (opts.blastops[l].compare("qcov") == 0)
					{
						// output % query coverage
						ss << "\t";
						ss.precision(3);
						ss << std::get<4>(miss_gap_match) * 100;
					}
					else if (opts.blastops[l].compare("qstrand") == 0)
					{
						// output strand
						ss << "\t";
						ss << strandmark;
					}
				}
				ss << std::endl;
			}//~blast tabular m8
		}
	} // ~iterate all alignments
	if (is_zip) {
		auto ret = vzlib_out[id].defstr(ss.str(), fsv[id]); // Z_STREAM_END | Z_OK - ok
		if (ret < Z_OK || ret > Z_STREAM_END) {
			ERR("Failed deflating readstring: ", ss.str(), " zlib status: ", ret);
		}
	}
	else {
		fsv[id] << ss.str();
	}
} // ~ ReportBlast::append