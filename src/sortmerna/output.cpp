/**
* @FILE: output.cpp
* @Created: Nov 26, 2017 Sun
* @brief Object for outputting results in various formats
* @parblock
* SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
* @copyright 2012-16 Bonsai Bioinformatics Research Group
* @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
*
* SortMeRNA is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* SortMeRNA is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with SortMeRNA.  If not, see <http://www.gnu.org/licenses/>.
* @endparblock
*
* @contributors Jenya Kopylova, jenya.kopylov@gmail.com
*               Laurent Noé, laurent.noe@lifl.fr
*               Pierre Pericard, pierre.pericard@lifl.fr
*               Daniel McDonald, wasade@gmail.com
*               Mikaël Salson, mikael.salson@lifl.fr
*               Hélène Touzet, helene.touzet@lifl.fr
*               Rob Knight, robknight@ucsd.edu
*/
#include "unistd.h"
#include <iomanip>
#include <fstream>
#include <cmath> // log, exp

#include "output.hpp"
#include "ThreadPool.hpp"
#include "kvdb.hpp"
#include "readsqueue.hpp"
#include "index.hpp"
#include "references.hpp"
#include "reader.hpp"
#include "processor.hpp"
#include "readstats.hpp"
#include "read.hpp"
#include "options.hpp"
#include "refstats.hpp"


// forward
void reportsJob(std::vector<Read> & reads, Runopts & opts, References & refs, Refstats & refstats, Output & output); // callback

void Output::init(Runopts & opts, Readstats & readstats)
{
	// attach pid to output files
	std::stringstream pidStr;
	if (opts.pid)
	{
		pidStr << getpid();
	}

	// associate the streams with reference sequence file names
	if (opts.filetype_ar.size() != 0)
	{
		if (opts.fastxout)
		{
			// fasta/fastq output
			fastaOutFile = opts.filetype_ar;
			if (opts.pid)
			{
				fastaOutFile.append("_");
				fastaOutFile.append(pidStr.str());
			}
			fastaOutFile.append(".");
			fastaOutFile.append(readstats.suffix);

			fastaout.open(fastaOutFile);
			fastaout.close();
		}

		if (opts.samout)
		{
			// sam output
			samoutFile = opts.filetype_ar;
			if (opts.pid)
			{
				samoutFile.append("_");
				samoutFile.append(pidStr.str());
			}
			samoutFile.append(".sam");
			samout.open(samoutFile);
			samout.close();
		}

		if (opts.blastout)
		{
			// blast output
			blastoutFile = opts.filetype_ar;
			if (opts.pid)
			{
				blastoutFile.append("_");
				blastoutFile.append(pidStr.str());
			}
			blastoutFile.append(".blast");
			blastout.open(blastoutFile);
			blastout.close();
		}

		// don't touch the log if only reports are generated
		if (opts.doLog && opts.alirep != Runopts::ALIGN_REPORT::report)
		{
			// statistics file output
			logfile = opts.filetype_ar;
			if (opts.pid)
			{
				logfile.append("_");
				logfile.append(pidStr.str());
			}
			logfile.append(".log");

			logstream.open(logfile);
			logstream.close();
		}

		if (opts.otumapout)
		{
			// OTU map output file
			std::ofstream otumap;
			otumapFile = opts.filetype_ar;
			if (opts.pid)
			{
				otumapFile.append("_");
				otumapFile.append(pidStr.str());
			}
			otumapFile.append("_otus.txt");
			otumap.open(otumapFile);
			otumap.close();
		}

		if (opts.de_novo_otu)
		{
			std::ofstream denovo_otu;
			denovo_otus_file = opts.filetype_ar;
			if (opts.pid)
			{
				denovo_otus_file.append("_");
				denovo_otus_file.append(pidStr.str());
			}
			denovo_otus_file.append("_denovo.");
			denovo_otus_file.append(readstats.suffix);

			denovo_otu.open(denovo_otus_file);
			denovo_otu.close();
		}
	}//~if ( ptr_filetype_ar != NULL ) 

	if (opts.filetype_or.size() != 0)
	{
		if (opts.fastxout)
		{
			// output stream for other reads
			std::ofstream fastaNonAlignOut;
			// add suffix database name to accepted reads file
			if (opts.pid)
			{
				opts.filetype_or += "_";
				opts.filetype_or += pidStr.str();
			}
			opts.filetype_or += ".";
			opts.filetype_or += readstats.suffix;
			// create the other reads file
			fastaNonAlignOut.open(opts.filetype_or);
			fastaNonAlignOut.close();
		}
	}
} // ~Output::init

void Output::report_blast
(
	Runopts & opts,
	Refstats & refstats,
	References & refs,
	Read & read
)
{
	const char MATCH = '|';
	const char MISMATCH = '*';
	const char INDEL = '-';

	uint32_t id = 0;
	uint32_t mismatches = 0;
	uint32_t gaps = 0;
	char strandmark = '+';
	bool flip03 = false; // flag to flip back to 03 on return

	if (read.is03) {
		read.flip34();
		flip03 = true;
	}

	// TODO: iterating all alignments for each reference part is an overhead. Alignments are pre-ordered, 
	//       so each new part corresponds to an index range of alignment vector. It's enough to loop 
	//       only that range.
	// iterate all alignments of the read
	for (int i = 0; i < read.hits_align_info.alignv.size(); ++i)
	{
		if (read.hits_align_info.alignv[i].index_num == refs.num 
			&& read.hits_align_info.alignv[i].part == refs.part)
		{
			uint32_t bitscore = (uint32_t)((float)((refstats.gumbel[refs.num].first)
				* (read.hits_align_info.alignv[i].score1) - std::log(refstats.gumbel[refs.num].second)) / (float)std::log(2));

			double evalue_score = (double)refstats.gumbel[refs.num].second
				* refstats.full_ref[refs.num]
				* refstats.full_read[refs.num]
				* std::exp(-refstats.gumbel[refs.num].first * read.hits_align_info.alignv[i].score1);

			std::string refseq = refs.buffer[read.hits_align_info.alignv[i].ref_seq].sequence;
			std::string ref_id = refs.buffer[read.hits_align_info.alignv[i].ref_seq].getId();

			if (read.hits_align_info.alignv[i].strand)
				strandmark = '+';
			else
				strandmark = '-';

			if (read.hits_align_info.alignv[i].strand == read.reversed)
				read.revIntStr(); // reverse if necessary

			// Blast-like pairwise alignment (only for aligned reads)
			if (opts.blastFormat == BlastFormat::REGULAR)
			{
				blastout << "Sequence ID: ";
				blastout << ref_id; // print only start of the header till first space
				blastout << std::endl;

				blastout << "Query ID: ";
				blastout << read.getSeqId();
				blastout << std::endl;

				blastout << "Score: " << read.hits_align_info.alignv[i].score1 << " bits (" << bitscore << ")\t";
				blastout.precision(3);
				blastout << "Expect: " << evalue_score << "\t";

				blastout << "strand: " << strandmark << std::endl << std::endl;

				if (read.hits_align_info.alignv[i].cigar.size() > 0)
				{
					uint32_t j, c = 0, left = 0, e = 0,
						qb = read.hits_align_info.alignv[i].ref_begin1,
						pb = read.hits_align_info.alignv[i].read_begin1;

					while (e < read.hits_align_info.alignv[i].cigar.size() || left > 0)
					{
						int32_t count = 0;
						int32_t q = qb;
						int32_t p = pb;
						blastout << "Target: ";
						blastout.width(8);
						blastout << q + 1 << "    ";
						// process CIGAR
						for (c = e; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
						{
							// 4 Low bits encode a Letter: M | D | S
							uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
							// 28 High bits encode the number of occurencies e.g. 34
							uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 1) blastout << INDEL; // mark indel
								else
								{
									blastout << nt_map[(int)refseq[q]];
									++q;
								}
								++count;
								if (count == 60) goto step2;
							}
						}
					step2:
						blastout << "    " << q << "\n";
						blastout.width(20);
						blastout << " ";
						q = qb;
						count = 0;
						for (c = e; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
						{
							//uint32_t letter = 0xf & *(a->cigar + c);
							uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 0)
								{
									if ((char)nt_map[(int)refseq[q]] == (char)nt_map[(int)read.isequence[p]]) blastout << MATCH; // mark match
									else blastout << MISMATCH; // mark mismatch
									++q;
									++p;
								}
								else
								{
									blastout << " ";
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
						blastout << "\nQuery: ";
						blastout.width(9);
						blastout << p + 1 << "    ";
						count = 0;
						for (c = e; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 2) blastout << INDEL; // mark indel
								else
								{
									blastout << nt_map[(int)read.isequence[p]];
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
						blastout << "    " << p << "\n\n";
					}
				}
			}
			// Blast tabular m8 + optional columns for CIGAR and query coverage
			else if (opts.blastFormat == BlastFormat::TABULAR)
			{
				// (1) Query ID
				blastout << read.getSeqId();

				// print null alignment for non-aligned read
				if (opts.print_all_reads && (read.hits_align_info.alignv.size() == 0))
				{
					blastout << "\t*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
					for (uint32_t l = 0; l < opts.blastops.size(); l++)
					{
						if (opts.blastops[l].compare("cigar") == 0)
							blastout << "\t*";
						else if (opts.blastops[l].compare("qcov") == 0)
							blastout << "\t0";
						else if (opts.blastops[l].compare("qstrand") == 0)
							blastout << "\t*";
						blastout << "\n";
					}
					return;
				}

				read.calcMismatchGapId(refs, i, mismatches, gaps, id);
				int32_t total_pos = mismatches + gaps + id;

				blastout << "\t";
				// (2) Subject
				blastout << ref_id << "\t";
				// (3) %id
				blastout.precision(3);
				blastout << (double)id / (mismatches + gaps + id) * 100 << "\t";
				// (4) alignment length
				blastout << (read.hits_align_info.alignv[i].read_end1 - read.hits_align_info.alignv[i].read_begin1 + 1) << "\t";
				// (5) mismatches
				blastout << mismatches << "\t";
				// (6) gap openings
				blastout << gaps << "\t";
				// (7) q.start
				blastout << read.hits_align_info.alignv[i].read_begin1 + 1 << "\t";
				// (8) q.end
				blastout << read.hits_align_info.alignv[i].read_end1 + 1 << "\t";
				// (9) s.start
				blastout << read.hits_align_info.alignv[i].ref_begin1 + 1 << "\t";
				// (10) s.end
				blastout << read.hits_align_info.alignv[i].ref_end1 + 1 << "\t";
				// (11) e-value
				blastout << evalue_score << "\t";
				// (12) bit score
				blastout << bitscore;
				// OPTIONAL columns
				for (uint32_t l = 0; l < opts.blastops.size(); l++)
				{
					// output CIGAR string
					if (opts.blastops[l].compare("cigar") == 0)
					{
						blastout << "\t";
						// masked region at beginning of alignment
						if (read.hits_align_info.alignv[i].read_begin1 != 0) blastout << read.hits_align_info.alignv[i].read_begin1 << "S";
						for (int c = 0; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
							blastout << length;
							if (letter == 0) blastout << "M";
							else if (letter == 1) blastout << "I";
							else blastout << "D";
						}

						uint32_t end_mask = read.sequence.length() - read.hits_align_info.alignv[i].read_end1 - 1;
						// output the masked region at end of alignment
						if (end_mask > 0) blastout << end_mask << "S";
					}
					// output % query coverage
					else if (opts.blastops[l].compare("qcov") == 0)
					{
						blastout << "\t";
						blastout.precision(3);
						double coverage = abs(read.hits_align_info.alignv[i].read_end1 - read.hits_align_info.alignv[i].read_begin1 + 1)
							/ read.hits_align_info.alignv[i].readlen;
						blastout << coverage * 100; // (double)align_len / readlen
					}
					// output strand
					else if (opts.blastops[l].compare("qstrand") == 0)
					{
						blastout << "\t";
						blastout << strandmark;
						//if (read.hits_align_info.alignv[i].strand) blastout << "+";
						//else blastout << "-";
					}
				}
				blastout << std::endl;
			}//~blast tabular m8
		}
	} // ~iterate all alignments

	if (flip03) read.flip34();
} // ~ Output::report_blast


void Output::writeSamHeader(Runopts & opts)
{
	samout << "@HD\tVN:1.0\tSO:unsorted\n";

	// TODO: this line is taken from "Index::load_stats". To be finished (20171215).
#if 0
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); index_num++)
	{
		//@SQ header
		if (opts.yes_SQ) acceptedsam << "@SQ\tSN:" << s << "\tLN:" << len_seq << "\n";
		// number of nucleotide sequences in the reference file
		uint32_t num_sq = 0;
		stats.read(reinterpret_cast<char*>(&num_sq), sizeof(uint32_t));

		// loop through each @SQ
		for (uint32_t j = 0; j < num_sq; j++)
		{
			// length of the sequence id
			uint32_t len_id = 0;
			stats.read(reinterpret_cast<char*>(&len_id), sizeof(uint32_t));
			// the sequence id string
			std::string s(len_id + 1, 0); // AK
			std::vector<char> vs(s.begin(), s.end());
			stats.read(reinterpret_cast<char*>(&vs[0]), sizeof(char)*len_id);
			// the length of the sequence itself
			uint32_t len_seq = 0;
			stats.read(reinterpret_cast<char*>(&len_seq), sizeof(uint32_t));
			//		 @SQ header
			if (opts.yes_SQ) acceptedsam << "@SQ\tSN:" << s << "\tLN:" << len_seq << "\n";
		} // ~for
	} // ~for
#endif
	samout << "@PG\tID:sortmerna\tVN:1.0\tCL:" << opts.cmdline << std::endl;

} // ~Output::writeSamHeader

void Output::report_sam
(
	Runopts & opts,
	References & refs,
	Read & read
)
{
	bool flip03 = false;
	if (read.is03) {
		read.flip34();
		flip03 = true;
	}

	//if (read.hits_align_info.alignv.size() == 0 && !opts.print_all_reads)
	//	return;

	// read did not align, output null string
	if (opts.print_all_reads && read.hits_align_info.alignv.size() == 0)
	{
		// (1) Query
		samout << read.getSeqId();
		samout << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		return;
	}

	// read aligned, output full alignment
	// iterate read alignments
	for (int i = 0; i < read.hits_align_info.alignv.size(); ++i)
	{
		if (read.hits_align_info.alignv[i].index_num == refs.num 
			&& read.hits_align_info.alignv[i].part == refs.part)
		{
			// (1) Query
			samout << read.getSeqId();
			// (2) flag Forward/Reversed
			if (!read.hits_align_info.alignv[i].strand) samout << "\t16\t";
			else samout << "\t0\t";
			// (3) Subject
			samout << refs.buffer[read.hits_align_info.alignv[i].ref_seq].getId();
			// (4) Ref start
			samout << "\t" << read.hits_align_info.alignv[i].ref_begin1 + 1;
			// (5) mapq
			samout << "\t" << 255 << "\t";
			// (6) CIGAR
			// output the masked region at beginning of alignment
			if (read.hits_align_info.alignv[i].read_begin1 != 0)
				samout << read.hits_align_info.alignv[i].read_begin1 << "S";

			for (int c = 0; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
			{
				uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
				uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
				samout << length;
				if (letter == 0) samout << "M";
				else if (letter == 1) samout << "I";
				else samout << "D";
			}

			uint32_t end_mask = read.sequence.size() - read.hits_align_info.alignv[i].read_end1 - 1;
			// output the masked region at end of alignment
			if (end_mask > 0) samout << end_mask << "S";
			// (7) RNEXT, (8) PNEXT, (9) TLEN
			samout << "\t*\t0\t0\t";
			// (10) SEQ

			if ( read.hits_align_info.alignv[i].strand == read.reversed ) // XNOR
				read.revIntStr();
			samout << read.get04alphaSeq();
			// (11) QUAL
			samout << "\t";
			// reverse-complement strand
			if (read.quality.size() > 0 && !read.hits_align_info.alignv[i].strand)
			{
				std::reverse(read.quality.begin(), read.quality.end());
				samout << read.quality;
			}
			else if (read.quality.size() > 0) // forward strand
			{
				samout << read.quality;
				// FASTA read
			}
			else samout << "*";

			// (12) OPTIONAL FIELD: SW alignment score generated by aligner
			samout << "\tAS:i:" << read.hits_align_info.alignv[i].score1;
			// (13) OPTIONAL FIELD: edit distance to the reference
			uint32_t mismatches = 0;
			uint32_t gaps = 0;
			uint32_t id = 0;
			read.calcMismatchGapId(refs, i, mismatches, gaps, id);
			samout << "\tNM:i:" << mismatches + gaps << "\n";
		}
	} // ~for read.alignments

	if (flip03) read.flip34();
} // ~Output::report_sam

/* 
 * prototype outputformats.cpp:report_fasta
 *
 * called on each Read or each 2 reads (if paired)
 *
 * @param reads: 1 or 2 reads (if paired)
 */
void Output::report_fasta(Runopts & opts, std::vector<Read> & reads)
{
	std::stringstream ss;

	// output accepted reads
	if (opts.fastxout && fastaout.is_open())
	{
		// pair-ended reads
		if (opts.pairedin || opts.pairedout)
		{
			// either both reads are accepted, or one is accepted and pairedin
			if ((reads[0].hit && reads[1].hit) ||
				((reads[0].hit || reads[1].hit) && opts.pairedin))
			{
				// output aligned read
				if (opts.fastxout)
				{
					for (Read read: reads)
					{
						fastaout << read.header << std::endl << read.sequence << std::endl;
						if (read.format == Format::FASTQ)
							fastaout << '+' << std::endl << read.quality << std::endl;
					}
				}
			}//~the read was accepted
		}//~if paired-in or paired-out
		else // regular or pair-ended reads don't need to go into the same file
		{
			// the read was accepted
			if (reads[0].hit)
			{
				// output aligned read
				if (opts.fastxout)
				{
					fastaout << reads[0].header << std::endl << reads[0].sequence << std::endl;
					if (reads[0].format == Format::FASTQ)
						fastaout << '+' << std::endl << reads[0].quality << std::endl;
				}
			} //~if read was accepted
		}//~if not paired-in or paired-out
	}//~if ( ptr_filetype_ar != NULL )

	// output other reads
	if (opts.fastxout && fastaNonAlignOut.is_open())
	{
		// pair-ended reads
		if (opts.pairedin || opts.pairedout)
		{
			// neither of the reads were accepted, or exactly one was accepted and pairedout_gv
			if ((!reads[0].hit && !reads[1].hit) ||
				((reads[0].hit ^ reads[1].hit) && opts.pairedout))
			{
				for (Read read : reads)
				{
					fastaNonAlignOut << read.header << std::endl << read.sequence << std::endl;
					if (read.format == Format::FASTQ)
						fastaNonAlignOut << '+' << std::endl << read.quality << std::endl;
				}
			}//~the read was accepted
		}//~if (pairedin_gv || pairedout_gv)
		else // output reads single
		{
			// the read was accepted
			if (!reads[0].hit)
			{
				fastaNonAlignOut << reads[0].header << std::endl << reads[0].sequence << std::endl;
				if (reads[0].format == Format::FASTQ)
					fastaNonAlignOut << '+' << std::endl << reads[0].quality << std::endl;
			}
		} // ~ if (!(pairedin_gv || pairedout_gv))
	} //~if ( opts.fastxout )  
} // ~Output::report_fasta

void Output::report_denovo(Runopts & opts, std::vector<Read> & reads)
{
	std::stringstream ss;

	// output reads with < id% alignment (passing E-value) for de novo clustering
	if (denovo_otus_file.size() != 0)
	{
		// pair-ended reads
		if (opts.pairedin || opts.pairedout)
		{
			// either both reads are accepted, or one is accepted and pairedin_gv
			if ( opts.pairedin && reads[0].hit && reads[1].hit && (reads[0].hit_denovo || reads[1].hit_denovo) )
			{
				// output aligned read
				for (Read read : reads)
					denovoreads << read.header << std::endl << read.sequence << std::endl;
			}//~the read was accepted
		}//~if paired-in or paired-out
		else // regular or pair-ended reads don't need to go into the same file
		{
			// the read was accepted
			if (reads[0].hit && reads[0].hit_denovo)
			{
				// output aligned read
				denovoreads << reads[0].header << std::endl << reads[0].sequence << std::endl;
			} //~if read was accepted
		}//~if not paired-in or paired-out
	}//~if ( denovo_otus_file set )
} // ~Output::report_denovo

void Output::report_biom(){

	biomout.open(biomfile, std::ios::in);

	if (biomout.is_open())
	{
		biomout << "\"id:\"null,";
		biomout << "\"format\": \"Biological Observation Matrix 1.0.0\",";
		biomout << "\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\"";
		biomout << "\"type\": \"OTU table\",";
		biomout << "\"generated_by\": \"SortMeRNA v2.0\",";
		biomout << "\"date\": \"\",";
		biomout << "\"rows\":[";
		biomout << "\"matrix_type\": \"sparse\",";
		biomout << "\"matrix_element_type\": \"int\",";
		biomout << "\"shape\":";
		biomout << "\"data\":";

		biomout.close();
	}
} // ~Output::report_biom

/**
 * open streams for writing
 */
void Output::openfiles(Runopts & opts)
{
	std::stringstream ss;

	if (opts.blastout && !blastout.is_open()) {
		blastout.open(blastoutFile);
		if (!blastout.good())
		{
			ss << "  " << RED << "ERROR" << COLOFF << ": could not open BLAST output file for writing.\n";
			std::cerr << ss.str(); ss.str("");
			exit(EXIT_FAILURE);
		}
	}

	if (opts.samout && !samout.is_open()) {
		samout.open(samoutFile);
		if (!samout.good())
		{
			ss << "  " << RED  << "ERROR" << COLOFF  << ": could not open SAM output file for writing.\n";
			std::cerr << ss.str(); ss.str("");
			exit(EXIT_FAILURE);
		}
	}

	if (opts.fastxout && !fastaout.is_open()) {
		fastaout.open(fastaOutFile, std::ios::app | std::ios::binary);
		if (!fastaout.good())
		{
			ss << "  " << RED << "ERROR" << COLOFF << ": could not open FASTA/Q output file for writing.\n";
			std::cerr << ss.str(); ss.str("");
			exit(EXIT_FAILURE);
		}
	}

	if (opts.fastxout && opts.filetype_or.size() != 0 && !fastaNonAlignOut.is_open())
	{
		fastaNonAlignOut.open(opts.filetype_or, std::ios::app | std::ios::binary);
		if (!fastaNonAlignOut.good())
		{
			ss << "  " << RED << "ERROR" << COLOFF << ": could not open FASTA/Q Non-aligned output file for writing." << std::endl;
			std::cerr << ss.str(); ss.str("");
			exit(EXIT_FAILURE);
		}
	}

	if (denovo_otus_file.size() != 0 && !denovoreads.is_open())
	{
		denovoreads.open(denovo_otus_file, std::ios::app | std::ios::binary);
		if (!denovoreads.good())
		{
			ss << "  " << RED << "ERROR" << denovo_otus_file << ": file " << COLOFF 
				<< " (denovoreads) could not be opened for writing." << std::endl;
			std::cerr << ss.str(); ss.str("");
			exit(EXIT_FAILURE);
		}
	}

	if (opts.doLog && logfile.size() > 0 && !logstream.is_open())
	{
		logstream.open(logfile, std::ofstream::binary | std::ofstream::app);
		if (!logstream.good())
		{
			ss << "  " << RED << "ERROR" << logfile << ": file " << COLOFF
				<< " (logfile) could not be opened." << std::endl;
			std::cerr << ss.str(); ss.str("");
			exit(EXIT_FAILURE);
		}
	}
} // ~Output::openfiles

void Output::closefiles()
{
	if (blastout.is_open()) blastout.close();
	if (samout.is_open()) samout.close();
	if (fastaout.is_open()) fastaout.close();
	if (fastaNonAlignOut.is_open()) fastaNonAlignOut.close();
	if (denovoreads.is_open()) denovoreads.close();
}

// called from main. TODO: move into a class?
void generateReports(Runopts & opts, Readstats & readstats, Output & output)
{
	int N_READ_THREADS = opts.num_read_thread_rep;
	int N_PROC_THREADS = opts.num_proc_thread_rep;
	int loopCount = 0; // counter of total number of processing iterations. TODO: no need here?
	std::stringstream ss;

	ss << "\tgenerateReports called. Thread: " << std::this_thread::get_id() << std::endl;
	std::cout << ss.str(); ss.str("");

	ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS);
	KeyValueDatabase kvdb(opts.kvdbPath);
	bool indb = readstats.restoreFromDb(kvdb);

	if (indb) {
		ss << __FILE__ << ":" << __LINE__ << " Restored Readstats from DB: " << indb << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	ReadsQueue readQueue("read_queue", QUEUE_SIZE_MAX, N_READ_THREADS); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", QUEUE_SIZE_MAX, N_PROC_THREADS); // Not used for Reports
	Refstats refstats(opts, readstats);
	References refs;

	output.openfiles(opts);
	if (opts.samout) output.writeSamHeader(opts);

	// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); ++index_num)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
		{
			ss << "\tLoading reference " << index_num << " part " << idx_part+1 << "/" << refstats.num_index_parts[index_num] << "  ... ";
			std::cout << ss.str(); ss.str("");
			auto starts = std::chrono::high_resolution_clock::now();
			refs.load(index_num, idx_part, opts, refstats);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts; // ~20 sec Debug/Win
			ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << " sec]\n";
			std::cout << ss.str(); ss.str("");

			starts = std::chrono::high_resolution_clock::now(); // index processing starts

			for (int i = 0; i < N_READ_THREADS; ++i)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), opts, readQueue, kvdb, loopCount));
			}

			// add processor jobs
			for (int i = 0; i < N_PROC_THREADS; ++i)
			{
				tpool.addJob(ReportProcessor("report_proc_" + std::to_string(i), readQueue, opts, refs, output, refstats, reportsJob));
			}
			++loopCount;
			tpool.waitAll(); // wait till processing is done on one index part
			refs.clear();
			writeQueue.reset(N_PROC_THREADS);
			readQueue.reset(N_READ_THREADS);

			elapsed = std::chrono::high_resolution_clock::now() - starts; // index processing done
			ss << "   Done reference " << index_num << " Part: " << idx_part + 1
				<< " Time: " << std::setprecision(2) << std::fixed << elapsed.count() << " sec\n";
			std::cout << ss.str(); ss.str("");
			if (!opts.blastout && !opts.samout)	goto done1;
		} // ~for(idx_part)
	} // ~for(index_num)
done1:
	output.closefiles();
	std::cout << "\tDone generateReports\n";
} // ~generateReports