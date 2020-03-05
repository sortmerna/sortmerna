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

#include "output.hpp"
#include "ThreadPool.hpp"
#include "kvdb.hpp"
#include "readsqueue.hpp"
#include "index.hpp"
#include "references.hpp"
#include "read_control.hpp"
#include "processor.hpp"
#include "readstats.hpp"
#include "read.hpp"
#include "options.hpp"
#include "refstats.hpp"


// forward
void reportsJob(std::vector<Read> & reads, Runopts & opts, References & refs, Refstats & refstats, Output & output); // callback

Summary::Summary():
	is_de_novo_otu(false), 
	is_otumapout(false), 
	total_reads(0), 
	total_reads_denovo_clustering(0),
	total_reads_mapped(0),
	total_reads_mapped_cov(0),
	min_read_len(0),
	max_read_len(0),
	all_reads_len(0),
	total_otu(0)
{}

Output::Output(Runopts& opts, Readstats& readstats)
{
	init(opts, readstats);
}
Output::~Output() { closefiles(); }

void Output::init(Runopts & opts, Readstats & readstats)
{
	summary.pid_str = std::to_string(getpid());

	// init aligned output files
	if (opts.is_fastx)
	{
		// fasta/q output  WORKDIR/out/aligned.fastq
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + summary.pid_str;
		}
		sfx += "." + readstats.suffix;
		size_t nfiles = opts.is_out2 ? 2 : 1;
		aligned_os.resize(nfiles);
		aligned_f.resize(nfiles);
		for (size_t i = 0; i < aligned_os.size(); ++i) {
			std::string sfx2 = "";
			if (opts.is_out2) {
				sfx2 = i == 0 ? "_fwd" : "_rev";
			}

			// test the file
			aligned_f[i] = opts.aligned_pfx.string() + sfx2 + sfx;
			std::cout << STAMP << "Testing file: " << std::filesystem::absolute(std::filesystem::path(aligned_f[i])) << std::endl;
			aligned_os[i].open(aligned_f[i]);
			aligned_os[i].close();
			if (!aligned_os[i]) {
				std::stringstream ss;
				ss << STAMP << "Failed operating stream on file " << aligned_f[i];
				ERR(ss.str());
				exit(EXIT_FAILURE);
			}
		}
	}

	if (opts.is_other && opts.is_fastx)
	{
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + summary.pid_str;
		}
		sfx += "." + readstats.suffix;
		size_t nfiles = opts.is_out2 ? 2 : 1;
		other_os.resize(nfiles);
		other_f.resize(nfiles);
		for (size_t i = 0; i < other_os.size(); ++i) {
			std::string sfx2 = "";
			if (opts.is_out2) {
				sfx2 = i == 0 ? "_fwd" : "_rev";
			}
			// test the file
			other_f[i] = opts.other_pfx.string() + sfx2 + sfx;
			std::cout << STAMP << "Testing file: " << std::filesystem::absolute(std::filesystem::path(other_f[i])) << std::endl;
			other_os[i].open(other_f[i]);
			other_os[i].close();
			if (!other_os[i]) {
				std::stringstream ss;
				ss << STAMP << "Failed operating stream on file " << other_f[i];
				ERR(ss.str());
				exit(EXIT_FAILURE);
			}
		}
	}

	if (opts.is_sam)
	{
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + summary.pid_str;
		}
		sfx += ".sam";
		sam_f = opts.aligned_pfx.string() + sfx;
		std::cout << STAMP << "Testing file: " << std::filesystem::absolute(std::filesystem::path(sam_f)) << std::endl;
		sam_os.open(sam_f);
		sam_os.close();
	}

	if (opts.is_blast)
	{
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + summary.pid_str;
		}
		sfx += ".blast";
		blast_f = opts.aligned_pfx.string() + sfx;
		std::cout << STAMP << "Testing file: " << std::filesystem::absolute(std::filesystem::path(blast_f)) << std::endl;
		blast_os.open(blast_f);
		blast_os.close();
		if (!blast_os) {
			std::stringstream ss;
			ss << STAMP << "Failed operating stream on file " << blast_f;
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
	}

	if (opts.is_otu_map)
	{
		// OTU map output file  WORKDIR/out/aligned_otus.txt
		std::ofstream otumap;
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + summary.pid_str;
		}
		sfx += "_otus.txt";
		otumap_f = opts.aligned_pfx.string() + sfx;
		std::cout << STAMP << "Testing file: " << std::filesystem::absolute(std::filesystem::path(otumap_f)) << std::endl;
		otumap.open(otumap_f);
		otumap.close();
	}

	if (opts.is_de_novo_otu)
	{
		std::ofstream denovo_otu;
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + summary.pid_str;
		}
		sfx += "_denovo." + readstats.suffix;
		denovo_otus_f = opts.aligned_pfx.string() + sfx;
		std::cout << STAMP << "Testing file: " << std::filesystem::absolute(std::filesystem::path(denovo_otus_f)) << std::endl;
		denovo_otu.open(denovo_otus_f);
		denovo_otu.close();
	}

	// don't touch the log if only reports are generated
	if (opts.is_log && opts.alirep != Runopts::ALIGN_REPORT::report)
	{
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + summary.pid_str;
		}
		sfx += ".log";
		log_f = opts.aligned_pfx.string() + sfx;
		std::cout << STAMP << "Testing file: " << std::filesystem::absolute(std::filesystem::path(log_f)) << std::endl;
		log_os.open(log_f);
		log_os.close();
		if (!log_os) {
			std::stringstream ss;
			ss << STAMP << "Failed operating stream on file " << log_f;
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
	}
} // ~Output::init

/**
 * called on each read => keep stream handle between calls 
 */
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

	if (read.is03) read.flip34();

	// TODO: iterating all alignments for each reference part is an overhead. Alignments are pre-ordered, 
	//       so each new part corresponds to an index range of alignment vector. It's enough to loop 
	//       only that range.
	// iterate all alignments of the read
	for (int i = 0; i < read.hits_align_info.alignv.size(); ++i)
	{
		if (read.hits_align_info.alignv[i].index_num == refs.num 
			&& read.hits_align_info.alignv[i].part == refs.part)
		{
			// (λ*S - ln(K))/ln(2)
			uint32_t bitscore = (uint32_t)((float)((refstats.gumbel[refs.num].first)
				* (read.hits_align_info.alignv[i].score1) - std::log(refstats.gumbel[refs.num].second)) / (float)std::log(2));

			// E = Kmn*exp(-λS)
			double evalue_score = (double)refstats.gumbel[refs.num].second
				* refstats.full_ref[refs.num]
				* refstats.full_read[refs.num]
				* std::exp(-refstats.gumbel[refs.num].first * read.hits_align_info.alignv[i].score1);

			std::string refseq = refs.buffer[read.hits_align_info.alignv[i].ref_seq].sequence;
			std::string ref_id = refs.buffer[read.hits_align_info.alignv[i].ref_seq].id;

			if (read.hits_align_info.alignv[i].strand)
				strandmark = '+';
			else
				strandmark = '-';

			if (read.hits_align_info.alignv[i].strand == read.reversed) // XNOR
				read.revIntStr(); // reverse if necessary

			// Blast-like pairwise alignment (only for aligned reads)
			if (opts.blastFormat == BlastFormat::REGULAR)
			{
				blast_os << "Sequence ID: ";
				blast_os << ref_id; // print only start of the header till first space
				blast_os << std::endl;

				blast_os << "Query ID: ";
				blast_os << read.getSeqId();
				blast_os << std::endl;

				blast_os << "Score: " << read.hits_align_info.alignv[i].score1 << " bits (" << bitscore << ")\t";
				blast_os.precision(3);
				blast_os << "Expect: " << evalue_score << "\t";

				blast_os << "strand: " << strandmark << std::endl << std::endl;

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
						blast_os << "Target: ";
						blast_os.width(8);
						blast_os << q + 1 << "    ";
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
								if (letter == 1) blast_os << INDEL; // mark indel
								else
								{
									blast_os << nt_map[(int)refseq[q]];
									++q;
								}
								++count;
								if (count == 60) goto step2;
							}
						}
					step2:
						blast_os << "    " << q << "\n";
						blast_os.width(20);
						blast_os << " ";
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
									if ((char)nt_map[(int)refseq[q]] == (char)nt_map[(int)read.isequence[p]]) blast_os << MATCH; // mark match
									else blast_os << MISMATCH; // mark mismatch
									++q;
									++p;
								}
								else
								{
									blast_os << " ";
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
						blast_os << "\nQuery: ";
						blast_os.width(9);
						blast_os << p + 1 << "    ";
						count = 0;
						for (c = e; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
							uint32_t l = (count == 0 && left > 0) ? left : length;
							for (j = 0; j < l; ++j)
							{
								if (letter == 2) blast_os << INDEL; // mark indel
								else
								{
									blast_os << nt_map[(int)read.isequence[p]];
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
						blast_os << "    " << p << "\n\n";
					}
				}
			}
			// Blast tabular m8 + optional columns for CIGAR and query coverage
			else if (opts.blastFormat == BlastFormat::TABULAR)
			{
				// (1) Query ID
				blast_os << read.getSeqId();

				// print null alignment for non-aligned read
				if (opts.is_print_all_reads && (read.hits_align_info.alignv.size() == 0))
				{
					blast_os << "\t*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
					for (uint32_t l = 0; l < opts.blastops.size(); l++)
					{
						if (opts.blastops[l].compare("cigar") == 0)
							blast_os << "\t*";
						else if (opts.blastops[l].compare("qcov") == 0)
							blast_os << "\t0";
						else if (opts.blastops[l].compare("qstrand") == 0)
							blast_os << "\t*";
						blast_os << "\n";
					}
					return;
				}

				read.calcMismatchGapId(refs, i, mismatches, gaps, id);
				int32_t total_pos = mismatches + gaps + id;

				blast_os << "\t";
				// (2) Subject
				blast_os << ref_id << "\t";
				// (3) %id
				blast_os.precision(3);
				blast_os << (double)id / (mismatches + gaps + id) * 100 << "\t";
				// (4) alignment length
				blast_os << (read.hits_align_info.alignv[i].read_end1 - read.hits_align_info.alignv[i].read_begin1 + 1) << "\t";
				// (5) mismatches
				blast_os << mismatches << "\t";
				// (6) gap openings
				blast_os << gaps << "\t";
				// (7) q.start
				blast_os << read.hits_align_info.alignv[i].read_begin1 + 1 << "\t";
				// (8) q.end
				blast_os << read.hits_align_info.alignv[i].read_end1 + 1 << "\t";
				// (9) s.start
				blast_os << read.hits_align_info.alignv[i].ref_begin1 + 1 << "\t";
				// (10) s.end
				blast_os << read.hits_align_info.alignv[i].ref_end1 + 1 << "\t";
				// (11) e-value
				blast_os << evalue_score << "\t";
				// (12) bit score
				blast_os << bitscore;
				// OPTIONAL columns
				for (uint32_t l = 0; l < opts.blastops.size(); l++)
				{
					// output CIGAR string
					if (opts.blastops[l].compare("cigar") == 0)
					{
						blast_os << "\t";
						// masked region at beginning of alignment
						if (read.hits_align_info.alignv[i].read_begin1 != 0) blast_os << read.hits_align_info.alignv[i].read_begin1 << "S";
						for (int c = 0; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
						{
							uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
							uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
							blast_os << length;
							if (letter == 0) blast_os << "M";
							else if (letter == 1) blast_os << "I";
							else blast_os << "D";
						}

						auto end_mask = read.sequence.length() - read.hits_align_info.alignv[i].read_end1 - 1;
						// output the masked region at end of alignment
						if (end_mask > 0) blast_os << end_mask << "S";
					}
					// output % query coverage
					else if (opts.blastops[l].compare("qcov") == 0)
					{
						blast_os << "\t";
						blast_os.precision(3);
						double coverage = abs(read.hits_align_info.alignv[i].read_end1 - read.hits_align_info.alignv[i].read_begin1 + 1)
							/ read.hits_align_info.alignv[i].readlen;
						blast_os << coverage * 100; // (double)align_len / readlen
					}
					// output strand
					else if (opts.blastops[l].compare("qstrand") == 0)
					{
						blast_os << "\t";
						blast_os << strandmark;
						//if (read.hits_align_info.alignv[i].strand) blastout << "+";
						//else blastout << "-";
					}
				}
				blast_os << std::endl;
			}//~blast tabular m8
		}
	} // ~iterate all alignments
} // ~ Output::report_blast


void Output::writeSamHeader(Runopts & opts)
{
	sam_os << "@HD\tVN:1.0\tSO:unsorted\n";

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
	sam_os << "@PG\tID:sortmerna\tVN:1.0\tCL:" << opts.cmdline << std::endl;

} // ~Output::writeSamHeader

void Output::report_sam
(
	Runopts & opts,
	References & refs,
	Read & read
)
{
	if (read.is03) read.flip34();

	//if (read.hits_align_info.alignv.size() == 0 && !opts.print_all_reads)
	//	return;

	// read did not align, output null string
	if (opts.is_print_all_reads && read.hits_align_info.alignv.size() == 0)
	{
		// (1) Query
		sam_os << read.getSeqId();
		sam_os << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
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
			sam_os << read.getSeqId();
			// (2) flag Forward/Reversed
			if (!read.hits_align_info.alignv[i].strand) sam_os << "\t16\t";
			else sam_os << "\t0\t";
			// (3) Subject
			sam_os << refs.buffer[read.hits_align_info.alignv[i].ref_seq].id;
			// (4) Ref start
			sam_os << "\t" << read.hits_align_info.alignv[i].ref_begin1 + 1;
			// (5) mapq
			sam_os << "\t" << 255 << "\t";
			// (6) CIGAR
			// output the masked region at beginning of alignment
			if (read.hits_align_info.alignv[i].read_begin1 != 0)
				sam_os << read.hits_align_info.alignv[i].read_begin1 << "S";

			for (int c = 0; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
			{
				uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
				uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
				sam_os << length;
				if (letter == 0) sam_os << "M";
				else if (letter == 1) sam_os << "I";
				else sam_os << "D";
			}

			auto end_mask = read.sequence.size() - read.hits_align_info.alignv[i].read_end1 - 1;
			// output the masked region at end of alignment
			if (end_mask > 0) sam_os << end_mask << "S";
			// (7) RNEXT, (8) PNEXT, (9) TLEN
			sam_os << "\t*\t0\t0\t";
			// (10) SEQ

			if ( read.hits_align_info.alignv[i].strand == read.reversed ) // XNOR
				read.revIntStr();
			sam_os << read.get04alphaSeq();
			// (11) QUAL
			sam_os << "\t";
			// reverse-complement strand
			if (read.quality.size() > 0 && !read.hits_align_info.alignv[i].strand)
			{
				std::reverse(read.quality.begin(), read.quality.end());
				sam_os << read.quality;
			}
			else if (read.quality.size() > 0) // forward strand
			{
				sam_os << read.quality;
				// FASTA read
			}
			else sam_os << "*";

			// (12) OPTIONAL FIELD: SW alignment score generated by aligner
			sam_os << "\tAS:i:" << read.hits_align_info.alignv[i].score1;
			// (13) OPTIONAL FIELD: edit distance to the reference
			uint32_t mismatches = 0;
			uint32_t gaps = 0;
			uint32_t id = 0;
			read.calcMismatchGapId(refs, i, mismatches, gaps, id);
			sam_os << "\tNM:i:" << mismatches + gaps << "\n";
		}
	} // ~for read.alignments
} // ~Output::report_sam

/* 
 * called on each Read or each 2 reads (if paired)
 * writes both aliged.fasta and other.fasta (if is_other is selected)
 *
 * @param reads: 1 or 2 (paired) reads
 */
void Output::report_fasta(Runopts & opts, std::vector<Read> & reads)
{
	std::stringstream ss;

	// output accepted reads
	if (opts.is_fastx)
	{
		//bool is_paired = opts.readfiles.size() == 2; // paired reads
		if (opts.is_paired)
		{
			if (opts.is_paired_in) {
				// if Either is aligned -> aligned
				if (reads[0].is_hit || reads[1].is_hit) {
					// validate the reads are paired in case of two reads files
					if (opts.readfiles.size() == 2 
						&& (reads[0].read_num != reads[1].read_num 
							|| reads[0].readfile_idx == reads[1].readfile_idx)) 
					{
						ss << STAMP << "Paired validation failed: reads[0].read_num = " << reads[0].read_num
							<< " reads[1].read_num = " << reads[0].read_num
							<< " reads[0].readfile_idx = " << reads[0].readfile_idx
							<< " reads[1].readfile_idx = " << reads[1].readfile_idx;
						ERR(ss.str());
						exit(EXIT_FAILURE);
					}
					
					// reads[0]
					for (size_t i = 0; i < reads.size(); ++i)
					{
						if (opts.is_out2) {
							write_a_read(aligned_os[i], reads[i]); // fwd and rev go into different files
						}
						else {
							write_a_read(aligned_os[0], reads[i]); // fwd and rev go into the same file
						}
					}
				}
				else if (opts.is_other) {
					for (size_t i = 0; i < reads.size(); ++i)
					{
						if (opts.is_out2) {
							write_a_read(other_os[i], reads[i]); // fwd and rev go into different files
						}
						else {
							write_a_read(other_os[0], reads[i]); // fwd and rev go into the same file
						}
					}
				}
			}
			else if (opts.is_paired_out) {
				// if Both aligned -> aligned
				if (reads[0].is_hit && reads[1].is_hit) {
					for (size_t i = 0; i < reads.size(); ++i)
					{
						if (opts.is_out2) {
							write_a_read(aligned_os[i], reads[i]); // fwd and rev go into different files
						}
						else {
							write_a_read(aligned_os[0], reads[i]); // fwd and rev go into the same file
						}
					}
				}
				// both non-aligned and is_other -> other
				else if (opts.is_other) {
					for (size_t i = 0; i < reads.size(); ++i)
					{
						if (opts.is_out2) {
							write_a_read(other_os[i], reads[i]); // fwd and rev go into different files
						}
						else {
							write_a_read(other_os[0], reads[i]); // fwd and rev go into the same file
						}
					}
				}
			}
			else {
				// Neither 'paired_in' nor 'paired_out' specified -> only aligned
				// reads go to aligned file, and only non-aligned to other file
				for (size_t i = 0; i < reads.size(); ++i)
				{
					if (reads[i].is_hit) {
						if (opts.is_out2) {
							write_a_read(aligned_os[i], reads[i]); // fwd and rev go into different files
						}
						else {
							write_a_read(aligned_os[0], reads[i]); // fwd and rev go into the same file
						}
					}
					else if (opts.is_other) {
						if (opts.is_out2) {
							write_a_read(other_os[i], reads[i]);
						}
						else {
							write_a_read(other_os[0], reads[i]); // fwd and rev go into the same file
						}
					}
				}
			}
		}//~if paired
		else // non-paired
		{
			// the read was accepted - output
			if (reads[0].is_hit)
			{
				write_a_read(aligned_os[0], reads[0]);
			} //~if read was accepted
			else if (opts.is_other) {
				write_a_read(other_os[0], reads[0]);
			}
		}//~if not paired-in or paired-out
	}//~if is_fastx 
} // ~Output::report_fasta

void Output::report_denovo(Runopts & opts, std::vector<Read> & reads)
{
	std::stringstream ss;

	// output reads with < id% alignment (passing E-value) for de novo clustering
	if (denovo_otus_f.size() != 0)
	{
		// pair-ended reads
		if (opts.is_paired_in || opts.is_paired_out)
		{
			// either both reads are accepted, or one is accepted and pairedin_gv
			if ( opts.is_paired_in && reads[0].is_hit && reads[1].is_hit && (reads[0].is_denovo || reads[1].is_denovo) )
			{
				// output aligned read
				for (Read read : reads)
					denovo_os << read.header << std::endl << read.sequence << std::endl;
			}//~the read was accepted
		}//~if paired-in or paired-out
		else // regular or pair-ended reads don't need to go into the same file
		{
			// the read was accepted
			if (reads[0].is_hit && reads[0].is_denovo)
			{
				// output aligned read
				denovo_os << reads[0].header << std::endl << reads[0].sequence << std::endl;
			} //~if read was accepted
		}//~if not paired-in or paired-out
	}//~if ( denovo_otus_file set )
} // ~Output::report_denovo

void Output::report_biom(){

	biom_os.open(biom_f, std::ios::in);

	if (biom_os.is_open())
	{
		biom_os << "\"id:\"null,";
		biom_os << "\"format\": \"Biological Observation Matrix 1.0.0\",";
		biom_os << "\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\"";
		biom_os << "\"type\": \"OTU table\",";
		biom_os << "\"generated_by\": \"SortMeRNA v2.0\",";
		biom_os << "\"date\": \"\",";
		biom_os << "\"rows\":[";
		biom_os << "\"matrix_type\": \"sparse\",";
		biom_os << "\"matrix_element_type\": \"int\",";
		biom_os << "\"shape\":";
		biom_os << "\"data\":";

		biom_os.close();
	}
} // ~Output::report_biom

/**
 * open streams for writing
 */
void Output::openfiles(Runopts & opts)
{
	std::stringstream ss;

	if (opts.is_blast && !blast_os.is_open()) {
		blast_os.open(blast_f);
		if (!blast_os.good())
		{
			ss.str("");
			ss << STAMP << "Could not open BLAST output file: [" << blast_f << "] for writing.";
			ERR(ss.str()); 
			exit(EXIT_FAILURE);
		}
	}

	if (opts.is_sam && !sam_os.is_open()) {
		sam_os.open(sam_f);
		if (!sam_os.good())
		{
			ss.str("");
			ss << STAMP  << "Could not open SAM output file ["<< sam_f << "] for writing.";
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
	}

	if (opts.is_fastx) {
		for (size_t i = 0; i < aligned_os.size(); ++i) {
			if (!aligned_os[i].is_open()) {
				aligned_os[i].open(aligned_f[i], std::ios::app | std::ios::binary);
			}
			if (!aligned_os[i].good())
			{
				ss.str("");
				ss << STAMP << "Could not open FASTA/Q output file [" << aligned_f[i] << "] for writing.";
				ERR(ss.str());
				exit(EXIT_FAILURE);
			}
		}
	}

	if (opts.is_fastx && opts.is_other)
	{
		for (size_t i = 0; i < other_os.size(); ++i) {
			if (!other_os[i].is_open()) {
				other_os[i].open(other_f[i], std::ios::app | std::ios::binary);
			}
			if (!other_os[i].good())
			{
				ss.str("");
				ss << STAMP << "Could not open FASTA/Q Non-aligned output file [" << other_f[i] << "] for writing.";
				ERR(ss.str());
				exit(EXIT_FAILURE);
			}
		}
	}

	if (denovo_otus_f.size() != 0 && !denovo_os.is_open())
	{
		denovo_os.open(denovo_otus_f, std::ios::app | std::ios::binary);
		if (!denovo_os.good())
		{
			ss.str("");
			ss << STAMP  << "Could not open denovo otus: [" << denovo_otus_f << "] for writing.";
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
	}
} // ~Output::openfiles

void Output::closefiles()
{
	if (blast_os.is_open()) { blast_os.flush(); blast_os.close(); }
	if (sam_os.is_open()) { sam_os.flush(); sam_os.close(); }
	for (size_t i = 0; i < aligned_os.size(); ++i) {
		if (aligned_os[i].is_open()) { aligned_os[i].flush(); aligned_os[i].close(); }
	}
	for (size_t i = 0; i < other_os.size(); ++i) {
		if (other_os[i].is_open()) { other_os[i].flush(); other_os[i].close(); }
	}
	if (denovo_os.is_open()) { denovo_os.flush(); denovo_os.close(); }

	std::cout << STAMP << "Flushed and closed" << std::endl;
}

/** 
 * called from postProcess 
 */
void Output::writeLog(Runopts &opts, Refstats &refstats, Readstats &readstats)
{
	if (!log_os.is_open())
	{
		log_os.open(log_f, std::ofstream::binary | std::ofstream::app);
	}

	std::cout << STAMP << "Using Log file: " << std::filesystem::absolute(log_f) << std::endl;

	summary.cmd = opts.cmdline;
	summary.total_reads = readstats.all_reads_count;
	if (opts.is_de_novo_otu) {
		summary.is_de_novo_otu = opts.is_de_novo_otu;
		summary.total_reads_denovo_clustering = readstats.total_reads_denovo_clustering;
	}
	summary.total_reads_mapped = readstats.total_reads_aligned.load();
	summary.min_read_len = readstats.min_read_len.load();
	summary.max_read_len = readstats.max_read_len.load();
	summary.all_reads_len = readstats.all_reads_len;

	// stats by database
	for (uint32_t index_num = 0; index_num < opts.indexfiles.size(); index_num++)
	{
		auto pcn = (float)((float)readstats.reads_matched_per_db[index_num] / (float)readstats.all_reads_count) * 100;
		summary.db_matches.emplace_back(std::make_pair(opts.indexfiles[index_num].first, pcn));
	}

	if (opts.is_otu_map) {
		summary.is_otumapout = opts.is_otu_map;
		summary.total_reads_mapped_cov = readstats.total_reads_mapped_cov.load();
		summary.total_otu = readstats.otu_map.size();
	}

	std::stringstream ss;
	time_t q = time(0);
	struct tm *now = localtime(&q);
	ss << asctime(now);
	summary.timestamp = ss.str();

	log_os << summary.to_string(opts, refstats);
#if 0
	logstream << " Command: [" << opts.cmdline << "]\n\n";

	// output total number of reads
	logstream << " Results:\n";
	logstream << "    Total reads = " << readstats.all_reads_count << std::endl;
	if (opts.de_novo_otu)
	{
		// all reads that have read::hit_denovo == true
		logstream << "    Total reads for de novo clustering = " << readstats.total_reads_denovo_clustering << std::endl;
	}
	// output total non-rrna + rrna reads
	logstream << std::setprecision(2) << std::fixed
		<< "    Total reads passing E-value threshold = " << readstats.total_reads_mapped.load()
		<< " (" << (float)((float)readstats.total_reads_mapped.load() / (float)readstats.all_reads_count) * 100 << ")" << std::endl
		<< "    Total reads failing E-value threshold = "
		<< readstats.all_reads_count - readstats.total_reads_mapped.load()
		<< " (" << (1 - ((float)((float)readstats.total_reads_mapped.load() / (float)readstats.all_reads_count))) * 100 << ")" << std::endl
		<< "    Minimum read length = " << readstats.min_read_len.load() << std::endl
		<< "    Maximum read length = " << readstats.max_read_len.load() << std::endl
		<< "    Mean read length    = " << readstats.all_reads_len / readstats.all_reads_count << std::endl
		<< " By database:" << std::endl;

	// output stats by database
	for (uint32_t index_num = 0; index_num < opts.indexfiles.size(); index_num++)
	{
		auto pcn = (float)((float)readstats.reads_matched_per_db[index_num] / (float)readstats.all_reads_count) * 100;
		std::make_pair(opts.indexfiles[index_num].first, pcn);
		logstream << "    " << opts.indexfiles[index_num].first << "\t\t" << pcn << std::endl;
	}

	if (opts.otumapout)
	{
		logstream << " Total reads passing %%id and %%coverage thresholds = " << readstats.total_reads_mapped_cov.load() << std::endl;
		logstream << " Total OTUs = " << readstats.otu_map.size() << std::endl;
	}
	logstream << std::endl << " " << asctime(now) << std::endl;
#endif
	log_os.close();
} // ~Output::writeLog

void Output::write_a_read(std::ofstream& strm, Read& read)
{
	strm << read.header << std::endl << read.sequence << std::endl;
	if (read.format == Format::FASTQ)
		strm << '+' << std::endl << read.quality << std::endl;
}

std::string Summary::to_string(Runopts &opts, Refstats &refstats)
{
	std::stringstream ss;

	ss << " Command:\n    " << cmd << std::endl << std::endl;

	ss << " Process pid = " << pid_str << std::endl << std::endl;

	ss << " Parameters summary: " << std::endl;
	int idx = 0;
	for (auto ref : opts.indexfiles) {
		ss << "    Reference file: " << ref.first << std::endl;
		ss << "        Seed length = " << opts.seed_win_len << std::endl;
		ss << "        Pass 1 = " << opts.skiplengths[idx][0] 
				<< ", Pass 2 = " << opts.skiplengths[idx][1] 
				<< ", Pass 3 = " << opts.skiplengths[idx][2] << std::endl;
		ss << "        Gumbel lambda = " << refstats.gumbel[idx].first << std::endl;
		ss << "        Gumbel K = " << refstats.gumbel[idx].second << std::endl;
		ss << "        Minimal SW score based on E-value = " << refstats.minimal_score[idx] << std::endl;
		++idx;
	}
	ss << "    Number of seeds = " << opts.seed_hits << std::endl;
	ss << "    Edges = " << opts.edges << std::endl;
	ss << "    SW match = " << opts.match << std::endl;
	ss << "    SW mismatch = " << opts.mismatch << std::endl;
	ss << "    SW gap open penalty = " << opts.gap_open << std::endl;
	ss << "    SW gap extend penalty = " << opts.gap_extension << std::endl;
	ss << "    SW ambiguous nucleotide = " << opts.score_N << std::endl;
	ss << "    SQ tags are " << (opts.is_SQ ? "" : "not ") << "output" << std::endl;
	ss << "    Number of alignment processing threads = " << opts.num_proc_thread << std::endl;
	for (auto readf : opts.readfiles) {
		ss << "    Reads file: " << readf << std::endl;
	}
	ss << "    Total reads = " << total_reads << std::endl << std::endl;

	ss << " Results:" << std::endl;
	if (is_de_novo_otu)
	{
		// all reads that have read::hit_denovo == true
		ss << "    Total reads for de novo clustering = " << total_reads_denovo_clustering << std::endl;
	}
	// output total non-rrna + rrna reads
	ss << std::setprecision(2) << std::fixed
		<< "    Total reads passing E-value threshold = " << total_reads_mapped
		<< " (" << ((float)total_reads_mapped / (float)total_reads * 100) << ")" << std::endl
		<< "    Total reads failing E-value threshold = " << total_reads - total_reads_mapped
		<< " (" << (1 - ((float)((float)total_reads_mapped / (float)total_reads))) * 100 << ")" << std::endl
		<< "    Minimum read length = " << min_read_len << std::endl
		<< "    Maximum read length = " << max_read_len << std::endl
		<< "    Mean read length    = " << all_reads_len / total_reads << std::endl << std::endl;

	ss << " Coverage by database:" << std::endl;

	// output stats by database
	for (auto match: db_matches)
	{
		ss << "    " << match.first << "\t\t" << match.second << std::endl;
	}

	if (is_otumapout)
	{
		ss << " Total reads passing %%id and %%coverage thresholds = " << total_reads_mapped_cov << std::endl;
		ss << " Total OTUs = " << total_otu << std::endl;
	}

	ss << std::endl << " " << timestamp << std::endl;

	return ss.str();
} // ~Summary::to_string

// called from main. TODO: move into a class?
void generateReports(Runopts & opts, Readstats & readstats, Output & output, KeyValueDatabase &kvdb)
{
	int N_READ_THREADS = opts.num_read_thread_rep;
	int N_PROC_THREADS = opts.num_proc_thread_rep;
	std::stringstream ss;

	ss.str("");
	ss << "\n" << STAMP << "=== Report generation starts. Thread: " << std::this_thread::get_id() << " ===\n\n";
	std::cout << ss.str();

	ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS);
	bool indb = readstats.restoreFromDb(kvdb);

	if (indb) {
		ss.str("");
		ss << STAMP << "Restored Readstats from DB: " << indb << std::endl;
		std::cout << ss.str(); 
	}

	ReadsQueue readQueue("read_queue", opts.queue_size_max, N_READ_THREADS); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", opts.queue_size_max, N_PROC_THREADS); // Not used for Reports
	Refstats refstats(opts, readstats);
	References refs;

	output.openfiles(opts);
	if (opts.is_sam) output.writeSamHeader(opts);

	// loop through every reference file passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); ++index_num)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
		{
			ss << std::endl << STAMP << "Loading reference " 
				<< index_num << " part " << idx_part+1 << "/" << refstats.num_index_parts[index_num] << "  ... ";
			std::cout << ss.str(); ss.str("");

			auto starts = std::chrono::high_resolution_clock::now();

			refs.load(index_num, idx_part, opts, refstats);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts; // ~20 sec Debug/Win
			ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << " sec]" << std::endl;
			std::cout << ss.str(); ss.str("");

			starts = std::chrono::high_resolution_clock::now(); // index processing starts

			for (int i = 0; i < N_READ_THREADS; ++i)
			{
				tpool.addJob(ReadControl(opts, readQueue, kvdb));
			}

			// add processor jobs
			for (int i = 0; i < N_PROC_THREADS; ++i)
			{
				tpool.addJob(ReportProcessor("report_proc_" + std::to_string(i), readQueue, opts, refs, output, refstats, reportsJob));
			}
			tpool.waitAll(); // wait till processing is done on one index part
			refs.clear();
			writeQueue.reset(N_PROC_THREADS);
			readQueue.reset(N_READ_THREADS);

			elapsed = std::chrono::high_resolution_clock::now() - starts; // index processing done
			ss.str("");
			ss << STAMP << "Done reference " << index_num << " Part: " << idx_part + 1
				<< " Time: " << std::setprecision(2) << std::fixed << elapsed.count() << " sec" << std::endl;
			std::cout << ss.str();
			if (!opts.is_blast && !opts.is_sam)	break;;
		} // ~for(idx_part)
	} // ~for(index_num)

	std::cout << "\n" << STAMP << "=== Done Reports generation ===\n\n";
} // ~generateReports