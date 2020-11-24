/**
* @FILE: output.cpp
* @Created: Nov 26, 2017 Sun
* @brief Object for outputting results in various formats
* @parblock
* SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
* @copyright 2016-20 Clarity Genomics BVBA
* @copyright 2012-16 Bonsai Bioinformatics Research Group
* @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
* @copyright 2016-19 Clarity Genomics BVBA
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
#include <filesystem>
#include <functional> // std::ref

#include "options.hpp"
#include "output.hpp"
#include "references.hpp"
#include "readstats.hpp"
#include "processor.hpp"
#include "refstats.hpp"
#include "readsqueue.hpp"
#include "readfeed.hpp"

// forward
class Read;
struct Index;
class KeyValueDatabase;

Output::Output(Readfeed& readfeed, Runopts& opts, Readstats& readstats) : out_type(0)
{
	init(readfeed, opts, readstats);
}
Output::~Output() { closefiles(); }

void Output::init(Readfeed& readfeed, Runopts& opts, Readstats& readstats)
{
	calc_out_type(opts);
	set_num_out();

	if (opts.is_fastx)
		init_fastx(readfeed, opts);

	if (opts.is_blast)
		init_blast(readfeed, opts);

	if (opts.is_sam)
		init_sam(readfeed, opts);

	if (opts.is_denovo)
		init_denovo(readfeed, opts);
} // ~Output::init

/*
 * af.fq, ar.fq
*/
void Output::init_fastx(Readfeed& readfeed, Runopts& opts)
{
	auto pid_str = std::to_string(getpid());
	auto num_aln = opts.is_other ? num_out >> 1 : num_out; // 1 | 2 | 4
	auto num_aln_split = readfeed.num_splits * num_aln;
	auto num_oth_split = opts.is_other ? num_aln_split : 0;
	ofs_aligned.resize(num_aln_split);
	f_aligned.resize(num_aln_split);
	if (opts.is_other) {
		f_other.resize(num_oth_split);
		ofs_other.resize(num_oth_split);
	}
	// fasta/q output  WORKDIR/out/aligned_paired_fwd_0_PID.fq
	//                              pfx + sfx1 + sfx2 + sfx3 + sfx4 + ext
	for (int i = 0; i < readfeed.num_splits; ++i) {
		for (int j = 0, idx = 0, orig_idx = 0; j < num_aln; ++j) {
			std::string sfx1 = "";
			if (out_type == 0x44 || out_type == 0x46) { // apf, apr, asf, asr
				if (j == 0 || j == 1) sfx1 = "_paired";
				else if (j == 2 || j == 3) sfx1 = "_singleton";
			}
			//(out_type == 0x5C || out_type == 0x5E) { // af, ar, of, or
			std::string sfx2 = "";
			sfx2 = opts.is_out2 && orig_idx == 0 ? "_fwd" : "_rev";

			std::string sfx3 = "_" + std::to_string(i);
			std::string sfx4 = opts.is_pid ? "_" + pid_str : "";
			std::string orig_ext = readfeed.orig_files[orig_idx].isFastq ? ".fq" : ".fa";

			// test the file(s)
			idx = i * num_out + j;
			f_aligned[idx] = opts.aligned_pfx.string() + sfx1 + sfx2 + sfx3 + orig_ext; // e.g. aligned_paired_fwd_0_PID.fq
			INFO("Testing file: ", std::filesystem::absolute(std::filesystem::path(f_aligned[idx])));
			ofs_aligned[idx].open(f_aligned[idx]);
			ofs_aligned[idx].close();
			if (!ofs_aligned[idx]) {
				ERR("Failed operating stream on file ", f_aligned[idx]);
				exit(EXIT_FAILURE);
			}

			if (opts.is_other) {
				f_other[idx] = opts.other_pfx.string() + sfx1 + sfx2 + sfx3 + sfx4 + orig_ext; // e.g. other_paired_fwd_0_PID.fq
				INFO("Testing file: ", std::filesystem::absolute(std::filesystem::path(f_other[idx])));
				ofs_other[idx].open(f_other[idx]);
				ofs_other[idx].close();
				if (!ofs_other[idx]) {
					ERR("Failed operating stream on file ", f_other[idx]);
					exit(EXIT_FAILURE);
				}
			}
			orig_idx = readfeed.is_two_files ? orig_idx^=1: orig_idx; // flip
		}
	}
} // ~Output::init_fastx

void Output::init_blast(Readfeed& readfeed, Runopts& opts)
{
	f_blast.resize(readfeed.num_splits);
	ofs_blast.resize(readfeed.num_splits);
	// WORKDIR/out/aligned_0_PID.blast
	std::string ext = ".blast";
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + std::to_string(getpid()) : "";
		f_blast[i] = opts.aligned_pfx.string() + sfx1 + sfx2 + ext;
		INFO("Testing file: ", std::filesystem::absolute(std::filesystem::path(f_blast[i])));
		ofs_blast[i].open(f_blast[i]);
		ofs_blast[i].close();
		if (!ofs_blast[i]) {
			ERR("Failed stream on file ", f_blast[i]);
			exit(EXIT_FAILURE);
		}
	}
} // ~Output::init_blast

void Output::init_sam(Readfeed& readfeed, Runopts& opts)
{
	f_sam.resize(readfeed.num_splits);
	ofs_sam.resize(readfeed.num_splits);
	// WORKDIR/out/aligned_0_PID.sam
	std::string ext = ".sam";
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + std::to_string(getpid()) : "";
		f_sam[i] = opts.aligned_pfx.string() + sfx1 + sfx2 + ext;
		INFO("Testing file: ", std::filesystem::absolute(std::filesystem::path(f_sam[i])));
		ofs_sam[i].open(f_sam[i]);
		ofs_sam[i].close();
		if (!ofs_sam[i]) {
			ERR("Failed stream on file ", f_sam[i]);
			exit(EXIT_FAILURE);
		}
	}
}

void Output::init_denovo(Readfeed& readfeed, Runopts& opts)
{
	f_denovo.resize(readfeed.num_splits);
	ofs_denovo.resize(readfeed.num_splits);
	// WORKDIR/out/aligned_denovo_0_PID.fa
	std::string ext = ".sam";
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + std::to_string(getpid()) : "";
		f_denovo[i] = opts.aligned_pfx.string() + "_denovo" + sfx1 + sfx2 + ".fa";
		INFO("Testing file: ", std::filesystem::absolute(std::filesystem::path(f_denovo[i])));
		ofs_denovo[i].open(f_denovo[i]);
		ofs_denovo[i].close();
		if (!ofs_denovo[i]) {
			ERR("Failed stream on file ", f_denovo[i]);
			exit(EXIT_FAILURE);
		}
	}
}

/**
 * called on each read => keep stream handle between calls 
 */
void Output::report_blast(int id, Read& read, References& refs, Refstats& refstats, Runopts& opts)
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
				ofs_blast[id] << "Sequence ID: ";
				ofs_blast[id] << ref_id; // print only start of the header till first space
				ofs_blast[id] << std::endl;

				ofs_blast[id] << "Query ID: ";
				ofs_blast[id] << read.getSeqId();
				ofs_blast[id] << std::endl;

				ofs_blast[id] << "Score: " << read.alignment.alignv[i].score1 << " bits (" << bitscore << ")\t";
				ofs_blast[id].precision(3);
				ofs_blast[id] << "Expect: " << evalue_score << "\t";

				ofs_blast[id] << "strand: " << strandmark << std::endl << std::endl;

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
						ofs_blast[id] << "Target: ";
						ofs_blast[id].width(8);
						ofs_blast[id] << q + 1 << "    ";
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
								if (letter == 1) ofs_blast[id] << INDEL; // mark indel
								else
								{
									ofs_blast[id] << nt_map[(int)refseq[q]];
									++q;
								}
								++count;
								if (count == 60) goto step2;
							}
						}
					step2:
						ofs_blast[id] << "    " << q << "\n";
						ofs_blast[id].width(20);
						ofs_blast[id] << " ";
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
									if ((char)nt_map[(int)refseq[q]] == (char)nt_map[(int)read.isequence[p]]) ofs_blast[id] << MATCH; // mark match
									else ofs_blast[id] << MISMATCH; // mark mismatch
									++q;
									++p;
								}
								else
								{
									ofs_blast[id] << " ";
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
						ofs_blast[id] << "\nQuery: ";
						ofs_blast[id].width(9);
						ofs_blast[id] << p + 1 << "    ";
						count = 0;
						for (c = e; c < read.alignment.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 2) ofs_blast[id] << INDEL; // mark indel
								else
								{
									ofs_blast[id] << nt_map[(int)read.isequence[p]];
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
						ofs_blast[id] << "    " << p << "\n\n";
					}
				}
			}
			// Blast tabular m8 + optional columns for CIGAR and query coverage
			else if (opts.blastFormat == BlastFormat::TABULAR)
			{
				// (1) Query ID
				ofs_blast[id] << read.getSeqId();

				// print null alignment for non-aligned read
				if (opts.is_print_all_reads && (read.alignment.alignv.size() == 0))
				{
					ofs_blast[id] << "\t*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
					for (uint32_t l = 0; l < opts.blastops.size(); l++)
					{
						if (opts.blastops[l].compare("cigar") == 0)
							ofs_blast[id] << "\t*";
						else if (opts.blastops[l].compare("qcov") == 0)
							ofs_blast[id] << "\t0";
						else if (opts.blastops[l].compare("qstrand") == 0)
							ofs_blast[id] << "\t*";
						ofs_blast[id] << "\n";
					}
					return;
				}

				read.calcMismatchGapId(refs, i, mismatches, gaps, mid);
				int32_t total_pos = mismatches + gaps + mid;

				ofs_blast[id] << "\t";
				// (2) Subject
				ofs_blast[id] << ref_id << "\t";
				// (3) %id
				ofs_blast[id].precision(3);
				ofs_blast[id] << (double)mid / total_pos * 100 << "\t";
				// (4) alignment length
				ofs_blast[id] << (read.alignment.alignv[i].read_end1 - read.alignment.alignv[i].read_begin1 + 1) << "\t";
				// (5) mismatches
				ofs_blast[id] << mismatches << "\t";
				// (6) gap openings
				ofs_blast[id] << gaps << "\t";
				// (7) q.start
				ofs_blast[id] << read.alignment.alignv[i].read_begin1 + 1 << "\t";
				// (8) q.end
				ofs_blast[id] << read.alignment.alignv[i].read_end1 + 1 << "\t";
				// (9) s.start
				ofs_blast[id] << read.alignment.alignv[i].ref_begin1 + 1 << "\t";
				// (10) s.end
				ofs_blast[id] << read.alignment.alignv[i].ref_end1 + 1 << "\t";
				// (11) e-value
				ofs_blast[id] << evalue_score << "\t";
				// (12) bit score
				ofs_blast[id] << bitscore;
				// OPTIONAL columns
				for (uint32_t l = 0; l < opts.blastops.size(); l++)
				{
					// output CIGAR string
					if (opts.blastops[l].compare("cigar") == 0)
					{
						ofs_blast[id] << "\t";
						// masked region at beginning of alignment
						if (read.alignment.alignv[i].read_begin1 != 0) ofs_blast[id] << read.alignment.alignv[i].read_begin1 << "S";
						for (int c = 0; c < read.alignment.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
							ofs_blast[id] << length;
							if (letter == 0) ofs_blast[id] << "M";
							else if (letter == 1) ofs_blast[id] << "I";
							else ofs_blast[id] << "D";
						}

						auto end_mask = read.sequence.length() - read.alignment.alignv[i].read_end1 - 1;
						// output the masked region at end of alignment
						if (end_mask > 0) ofs_blast[id] << end_mask << "S";
					}
					// output % query coverage
					else if (opts.blastops[l].compare("qcov") == 0)
					{
						double coverage = (double)abs(read.alignment.alignv[i].read_end1 - read.alignment.alignv[i].read_begin1 + 1)
							/ read.alignment.alignv[i].readlen;

						ofs_blast[id] << "\t";
						ofs_blast[id].precision(3);
						ofs_blast[id] << coverage * 100; // (double)align_len / readlen
					}
					// output strand
					else if (opts.blastops[l].compare("qstrand") == 0)
					{
						ofs_blast[id] << "\t";
						ofs_blast[id] << strandmark;
						//if (read.alignment.alignv[i].strand) blastout << "+";
						//else blastout << "-";
					}
				}
				ofs_blast[id] << std::endl;
			}//~blast tabular m8
		}
	} // ~iterate all alignments
} // ~ Output::report_blast


void Output::writeSamHeader(Runopts& opts)
{
	ofs_sam[0] << "@HD\tVN:1.0\tSO:unsorted\n";

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
	ofs_sam[0] << "@PG\tID:sortmerna\tVN:1.0\tCL:" << opts.cmdline << std::endl;

} // ~Output::writeSamHeader

void Output::report_sam(int id, Read& read, References& refs, Runopts& opts)
{
	if (read.is03) read.flip34();

	// read did not align, output null string
	if (opts.is_print_all_reads && read.alignment.alignv.size() == 0)
	{
		// (1) Query
		ofs_sam[id] << read.getSeqId();
		ofs_sam[id] << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		return;
	}

	// read aligned, output full alignment
	// iterate read alignments
	for (int i = 0; i < read.alignment.alignv.size(); ++i)
	{
		if (read.alignment.alignv[i].index_num == refs.num 
			&& read.alignment.alignv[i].part == refs.part)
		{
			// (1) Query
			ofs_sam[id] << read.getSeqId();
			// (2) flag Forward/Reversed
			if (!read.alignment.alignv[i].strand) ofs_sam[id] << "\t16\t";
			else ofs_sam[id] << "\t0\t";
			// (3) Subject
			ofs_sam[id] << refs.buffer[read.alignment.alignv[i].ref_num].id;
			// (4) Ref start
			ofs_sam[id] << "\t" << read.alignment.alignv[i].ref_begin1 + 1;
			// (5) mapq
			ofs_sam[id] << "\t" << 255 << "\t";
			// (6) CIGAR
			// output the masked region at beginning of alignment
			if (read.alignment.alignv[i].read_begin1 != 0)
				ofs_sam[id] << read.alignment.alignv[i].read_begin1 << "S";

			for (int c = 0; c < read.alignment.alignv[i].cigar.size(); ++c)
			{
				uint32_t letter = 0xf & read.alignment.alignv[i].cigar[c];
				uint32_t length = (0xfffffff0 & read.alignment.alignv[i].cigar[c]) >> 4;
				ofs_sam[id] << length;
				if (letter == 0) ofs_sam[id] << "M";
				else if (letter == 1) ofs_sam[id] << "I";
				else ofs_sam[id] << "D";
			}

			auto end_mask = read.sequence.size() - read.alignment.alignv[i].read_end1 - 1;
			// output the masked region at end of alignment
			if (end_mask > 0) ofs_sam[id] << end_mask << "S";
			// (7) RNEXT, (8) PNEXT, (9) TLEN
			ofs_sam[id] << "\t*\t0\t0\t";
			// (10) SEQ

			if ( read.alignment.alignv[i].strand == read.reversed ) // XNOR
				read.revIntStr();
			ofs_sam[id] << read.get04alphaSeq();
			// (11) QUAL
			ofs_sam[id] << "\t";
			// reverse-complement strand
			if (read.quality.size() > 0 && !read.alignment.alignv[i].strand)
			{
				std::reverse(read.quality.begin(), read.quality.end());
				ofs_sam[id] << read.quality;
			}
			else if (read.quality.size() > 0) // forward strand
			{
				ofs_sam[id] << read.quality;
				// FASTA read
			}
			else ofs_sam[id] << "*";

			// (12) OPTIONAL FIELD: SW alignment score generated by aligner
			ofs_sam[id] << "\tAS:i:" << read.alignment.alignv[i].score1;
			// (13) OPTIONAL FIELD: edit distance to the reference
			uint32_t mismatches = 0;
			uint32_t gaps = 0;
			uint32_t mid = 0;
			read.calcMismatchGapId(refs, i, mismatches, gaps, mid);
			ofs_sam[id] << "\tNM:i:" << mismatches + gaps << "\n";
		}
	} // ~for read.alignments
} // ~Output::report_sam

/* 
 * called on each Read or each 2 reads (if paired)
 * writes both aliged.fasta and other.fasta (if is_other is selected)
 *
 * @param reads: 1 or 2 (paired) reads
 */
void Output::report_fastx(int id, std::vector<Read>& reads, Runopts& opts)
{
	if (reads.size() == 2)
	{
		if (opts.is_paired_in) {
			// if Either is aligned -> both to aligned
			if (reads[0].is_hit || reads[1].is_hit) {
				// validate the reads are paired in case of two reads files
				if (opts.readfiles.size() == 2 
					&& (reads[0].read_num != reads[1].read_num 
						|| reads[0].readfile_idx == reads[1].readfile_idx)) 
				{
					ERR("Paired validation failed: reads[0].id= ", reads[0].id, " reads[0].read_num = ", 
						reads[0].read_num, " reads[0].readfile_idx= ", reads[0].readfile_idx,
						" reads[1].id=", reads[1].id, " reads[1].read_num = ", reads[1].read_num, 
						" reads[1].readfile_idx = ", reads[1].readfile_idx);
					exit(EXIT_FAILURE);
				}
				
				// reads[0]
				for (int i = 0, idx = id*reads.size(); i < reads.size(); ++i)
				{
					if (opts.is_out2) {
						write_a_read(ofs_aligned[idx], reads[i]);
						++idx; // fwd and rev go into different files
					}
					else {
						write_a_read(ofs_aligned[idx], reads[i]); // fwd and rev go into the same file
					}
				}
			}
			else if (opts.is_other) {
				for (int i = 0, idx = id * reads.size(); i < reads.size(); ++i)
				{
					if (opts.is_out2) {
						write_a_read(ofs_other[idx], reads[i]); // fwd and rev go into different files
						++idx;
					}
					else {
						write_a_read(ofs_other[0], reads[i]); // fwd and rev go into the same file
					}
				}
			}
		}
		else if (opts.is_paired_out) {
			// if Both aligned -> aligned
			if (reads[0].is_hit && reads[1].is_hit) {
				for (int i = 0, idx = id * reads.size(); i < reads.size(); ++i)
				{
					if (opts.is_out2) {
						write_a_read(ofs_aligned[idx], reads[i]); // fwd and rev go into different files
						++idx;
					}
					else {
						write_a_read(ofs_aligned[idx], reads[i]); // fwd and rev go into the same file
					}
				}
			}
			// both non-aligned and is_other -> other
			else if (opts.is_other) {
				for (int i = 0, idx = id * reads.size(); i < reads.size(); ++i)
				{
					if (opts.is_out2) {
						write_a_read(ofs_other[idx], reads[i]);
						++idx; // fwd and rev go into different files
					}
					else {
						write_a_read(ofs_other[idx], reads[i]); // fwd and rev go into the same file
					}
				}
			}
		}
		else {
			// Neither 'paired_in' nor 'paired_out' specified:
			//   aligned     -> aligned file
			//   non-aligned -> other file
			auto idx = id * num_out;
			if (reads[0].is_hit && reads[1].is_hit) {
				if (opts.is_out2) {
					write_a_read(ofs_aligned[idx], reads[0]); // apf
					write_a_read(ofs_aligned[idx+1], reads[1]); // apr
				}
				else {
					write_a_read(ofs_aligned[idx], reads[0]); // ap
					write_a_read(ofs_aligned[idx], reads[1]); // ap
				}
			}
			else if (reads[0].is_hit) { // fwd hit
				if (opts.is_out2) {
					write_a_read(ofs_aligned[idx+2], reads[0]); // asf
					if (opts.is_other) {
						write_a_read(ofs_other[idx + 2], reads[1]); // osf
					}
				}
				else {
					write_a_read(ofs_aligned[idx], reads[0]); // as
					if (opts.is_other) {
						write_a_read(ofs_other[idx], reads[1]); // os
					}
				}
			}
			else if (reads[1].is_hit) { // rev hit
				if (opts.is_out2) {
					write_a_read(ofs_aligned[idx+3], reads[0]); // asr
					if (opts.is_other) {
						write_a_read(ofs_other[idx+3], reads[1]); // osr
					}
				}
				else {
					write_a_read(ofs_aligned[idx+1], reads[0]); // as
					if (opts.is_other) {
						write_a_read(ofs_other[idx+1], reads[1]); // os
					}
				}
			}
			else { // none is hit
				if (opts.is_other) {
					if (opts.is_out2) {
						write_a_read(ofs_other[idx], reads[0]); // opf
						write_a_read(ofs_other[idx + 1], reads[1]); // opr
					}
					else {
						write_a_read(ofs_other[idx], reads[0]); // op
						write_a_read(ofs_other[idx], reads[1]); // op
					}
				}
			}
		}
	}//~if paired
	// non-paired
	else
	{
		// the read was accepted - output
		if (reads[0].is_hit)
		{
			write_a_read(ofs_aligned[0], reads[0]);
		}
		else if (opts.is_other) {
			write_a_read(ofs_other[0], reads[0]);
		}
	}
} // ~Output::report_fasta

/* 
 * output the 'de novo' reads i.e. reads that pass SW & fail (%Cov & %ID)
 * called on each read or a pair
 */
void Output::report_denovo(int id, Read& read, Runopts& opts)
{
	if (read.is_denovo)
		ofs_denovo[id] << read.header << std::endl << read.sequence << std::endl;
} // ~Output::report_denovo

void Output::report_biom(){

	ofs_biom.open(f_biom, std::ios::in);

	if (ofs_biom.is_open())
	{
		ofs_biom << "\"id:\"null,";
		ofs_biom << "\"format\": \"Biological Observation Matrix 1.0.0\",";
		ofs_biom << "\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\"";
		ofs_biom << "\"type\": \"OTU table\",";
		ofs_biom << "\"generated_by\": \"SortMeRNA v2.0\",";
		ofs_biom << "\"date\": \"\",";
		ofs_biom << "\"rows\":[";
		ofs_biom << "\"matrix_type\": \"sparse\",";
		ofs_biom << "\"matrix_element_type\": \"int\",";
		ofs_biom << "\"shape\":";
		ofs_biom << "\"data\":";

		ofs_biom.close();
	}
} // ~Output::report_biom

/**
 * open streams for writing
 */
void Output::openfiles(Runopts& opts)
{
	if (opts.is_fastx) {
		for (size_t i = 0; i < f_aligned.size(); ++i) {
			// aligned
			if (!ofs_aligned[i].is_open()) {
				ofs_aligned[i].open(f_aligned[i], std::ios::app | std::ios::binary);
			}
			if (!ofs_aligned[i].good()) {
				ERR("Could not open FASTA/Q output file [", f_aligned[i], "] for writing.");
				exit(EXIT_FAILURE);
			}
			// other
			if (opts.is_other) {
				if (!ofs_other[i].is_open()) {
					ofs_other[i].open(f_other[i], std::ios::app | std::ios::binary);
				}
				if (!ofs_other[i].good())
				{
					ERR("Could not open FASTA/Q Non-aligned output file [", f_other[i], "] for writing.");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	if (opts.is_blast) {
		for (int i = 0; i < ofs_blast.size(); ++i) {
			if (!ofs_blast[i].is_open()) {
				ofs_blast[i].open(f_blast[i]);
				if (!ofs_blast[i].good())
				{
					ERR("Could not open BLAST output file: [", f_blast[i], "] for writing.");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	if (opts.is_sam) {
		for (int i = 0; i < ofs_sam.size(); ++i) {
			if (!ofs_sam[i].is_open()) {
				ofs_sam[i].open(f_sam[i]);
				if (!ofs_sam[i].good()) {
					ERR("Could not open SAM output file ", f_sam[i], " for writing.");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	if (opts.is_denovo)	{
		for (int i = 0; i < ofs_denovo.size(); ++i) {
			if (!ofs_denovo[i].is_open()) {
				ofs_denovo[i].open(f_denovo[i], std::ios::app | std::ios::binary);
				if (!ofs_denovo[i].good()) {
					ERR("Could not open: ", f_denovo[i], " for writing.");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
} // ~Output::openfiles

void Output::closefiles()
{
	// fastx
	for (int i = 0; i < ofs_aligned.size(); ++i) {
		if (ofs_aligned[i].is_open()) { ofs_aligned[i].flush(); ofs_aligned[i].close(); }
	}
	for (int i = 0; i < ofs_other.size(); ++i) {
		if (ofs_other[i].is_open()) { ofs_other[i].flush(); ofs_other[i].close(); }
	}
	// blast
	for (int i = 0; i < ofs_blast.size(); ++i) {
		if (ofs_blast[i].is_open()) { 
			ofs_blast[i].flush(); ofs_blast[i].close(); 
		}
	}
	// sam
	for (int i = 0; i < ofs_sam.size(); ++i) {
		if (ofs_sam[i].is_open()) {
			ofs_sam[i].flush(); ofs_sam[i].close();
		}
	}
	// denovo
	for (int i = 0; i < ofs_denovo.size(); ++i) {
		if (ofs_denovo[i].is_open()) { 
			ofs_denovo[i].flush(); ofs_denovo[i].close(); 
		}
	}

	INFO("Flushed and closed");
}

void Output::write_a_read(std::ofstream& strm, Read& read)
{
	strm << read.header << std::endl << read.sequence << std::endl;
	if (read.format == BIO_FORMAT::FASTQ)
		strm << '+' << std::endl << read.quality << std::endl;
}

void Output::calc_out_type(Runopts& opts)
{
	if (opts.readfiles.size() == 1) out_type |= mask_1_file;
	if (opts.is_paired) out_type |= mask_paired;
	if (opts.readfiles.size() == 2) out_type |= mask_2_file;
	if (opts.is_other) out_type |= mask_other;
	if (opts.is_paired_in) out_type |= mask_paired_in;
	if (opts.is_paired_out) out_type |= mask_paired_out;
	if (opts.is_out2) out_type |= mask_out2;

	INFO("Output type: ", out_type);
}

void Output::set_num_out()
{
	if (out_type == 0x04 || out_type == 0x06) num_out = 2; // ap, as
	else if (out_type == 0x14 || out_type == 0x16 
		|| out_type == 0x24 || out_type == 0x26) num_out = 1; // a
	else if (out_type == 0x44 || out_type == 0x46) num_out = 4; // apf, apr, asf, asr   'af.size != ar.size' in general. Only aligned reads.
	else if (out_type == 0x4C || out_type == 0x4E) num_out = 8;
	else if (out_type == 0x5C || out_type == 0x5E) num_out = 4; // af, ar, of, or
}

void Output::merge_fastx(int num_splits)
{
	for (int i = 0; i < num_out; ++i) {
		std::ofstream ofs(f_aligned[i], std::ios_base::app | std::ios_base::binary);
		if (!ofs.is_open()) {
			ERR("failed to open for writing: ", f_aligned[i]);
			exit(1);
		}
		for (int j = 1; j < num_splits; ++j) {
			auto idx = i + j*num_out;
			std::ifstream ifs(f_aligned[idx], std::ios_base::out | std::ios_base::binary);
			if (ifs.is_open()) {
				ofs << ifs.rdbuf();
				INFO("merged ", f_aligned[idx], " -> ", f_aligned[i]);
				ifs.close();
				std::filesystem::remove(f_aligned[idx]);
				INFO("deleted ", f_aligned[idx]);
			}
			else {
				ERR("failed to open for reading: ", f_aligned[idx]);
				exit(1);
			}
		}
		ofs.close();
	}
}

void Output::merge_blast(int num_splits)
{
	std::ofstream ofs(f_blast[0], std::ios_base::app | std::ios_base::binary);
	if (!ofs.is_open()) {
		ERR("failed to open for writing: ", f_blast[0]);
		exit(1);
	}
	for (int i = 1; i < num_splits; ++i) {
		std::ifstream ifs(f_blast[i], std::ios_base::out | std::ios_base::binary);
		if (ifs.is_open()) {
			ofs << ifs.rdbuf();
			INFO("merged ", f_blast[i], " -> ", f_blast[0]);
			ifs.close();
			std::filesystem::remove(f_blast[i]);
			INFO("deleted ", f_blast[i]);
		}
		else {
			ERR("failed to open for reading: ", f_blast[i]);
			exit(1);
		}
	}
}

void Output::merge_sam(int num_splits)
{
	std::ofstream ofs(f_sam[0], std::ios_base::app | std::ios_base::binary);
	if (!ofs.is_open()) {
		ERR("failed to open for writing: ", f_sam[0]);
		exit(1);
	}
	for (int i = 1; i < num_splits; ++i) {
		std::ifstream ifs(f_sam[i], std::ios_base::out | std::ios_base::binary);
		if (ifs.is_open()) {
			ofs << ifs.rdbuf();
			INFO("merged ", f_sam[i], " -> ", f_sam[0]);
			ifs.close();
			std::filesystem::remove(f_sam[i]);
			INFO("deleted ", f_sam[i]);
		}
		else {
			ERR("failed to open for reading: ", f_sam[i]);
			exit(1);
		}
	}
}

void Output::merge_denovo(int num_splits)
{
	std::ofstream ofs(f_denovo[0], std::ios_base::app | std::ios_base::binary);
	if (!ofs.is_open()) {
		ERR("failed to open for writing: ", f_denovo[0]);
		exit(1);
	}
	for (int i = 1; i < num_splits; ++i) {
		std::ifstream ifs(f_denovo[i], std::ios_base::out | std::ios_base::binary);
		if (ifs.is_open()) {
			ofs << ifs.rdbuf();
			INFO("merged ", f_denovo[i], " -> ", f_denovo[0]);
			ifs.close();
			std::filesystem::remove(f_denovo[i]);
			INFO("deleted ", f_denovo[i]);
		}
		else {
			ERR("failed to open for reading: ", f_denovo[i]);
			exit(1);
		}
	}
}

void report(int id,
	Readfeed& readfeed,
	References& refs,
	Refstats& refstats,
	KeyValueDatabase& kvdb,
	Output& output,
	Runopts& opts)
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
		auto idx = id * readfeed.num_orig_files;
		for (std::size_t i = 0; i < num_reads; ++i)
		{
			if (readfeed.next(idx, readstr))
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
			idx = i == 0 ? idx + 1 : idx - 1; // switch fwd-rev
		}

		if (!isDone) {
			if (reads.back().isEmpty || !reads.back().isValid) {
				++num_invalid;
			}
			else {
				// only needs one loop through all reads - reference file is not used
				if (opts.is_fastx && refs.num == 0 && refs.part == 0)
				{
					output.report_fastx(id, reads, opts);
				}

				for (int i = 0; i < reads.size(); ++i)
				{
					if (opts.is_blast) output.report_blast(id, reads[i], refs, refstats, opts);
					if (opts.is_sam) output.report_sam(id, reads[i], refs, opts);
					// only needs one loop through all reads, no reference file dependency
					if (opts.is_denovo && refs.num == 0 && refs.part == 0) output.report_denovo(id, reads[i], opts);
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
		readfeed.init_reading(); // prepare readfeed
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
			INFO("done. Elapsed sec [", elapsed.count(), "]");

			starts = std::chrono::high_resolution_clock::now(); // index processing starts

			// start processor
			if (opts.feed_type == FEED_TYPE::SPLIT_READS) {
				for (int i = 0; i < nthreads; ++i) {
					tpool.emplace_back(std::thread(report, i, std::ref(readfeed),
						std::ref(refs), std::ref(refstats), std::ref(kvdb), std::ref(output), std::ref(opts)));
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

			refs.unload();
			INFO_MEM("References unloaded.");
			tpool.clear();
			// rewind for the next index
			readfeed.rewind_in();
			readfeed.init_vzlib_in();

			if (!opts.is_blast && !opts.is_sam)	break;;
		} // ~for(idx_part)
	} // ~for(ref_idx)

	output.closefiles();
	if (opts.is_fastx)
		output.merge_fastx(readfeed.num_splits);
	if (opts.is_blast)
		output.merge_blast(readfeed.num_splits);
	if (opts.is_sam)
		output.merge_sam(readfeed.num_splits);
	if (opts.is_denovo)
		output.merge_denovo(readfeed.num_splits);

	INFO("=== Done Reports ===\n");
} // ~writeReports