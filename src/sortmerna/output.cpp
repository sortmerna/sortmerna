/**
* FILE: output.cpp
* Created: Nov 26, 2017 Sun
*/
#include "unistd.h"
#include <iomanip>

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
	char pidStr[4000];
	if (opts.pid)
	{
		int32_t pid = _getpid();
		sprintf(pidStr, "%d", pid);
	}

	// associate the streams with reference sequence file names
	if (opts.filetype_ar.size() != 0)
	{
		if (opts.fastxout)
		{
			// fasta/fastq output
			acceptedstrings = opts.filetype_ar;
			if (opts.pid)
			{
				acceptedstrings.append("_");
				acceptedstrings.append(pidStr);
			}
			acceptedstrings.append(".");
			acceptedstrings.append(readstats.suffix.c_str());

			acceptedreads.open(acceptedstrings);
			acceptedreads.close();
		}

		if (opts.samout)
		{
			// sam output
			acceptedstrings_sam = opts.filetype_ar;
			if (opts.pid)
			{
				acceptedstrings_sam.append("_");
				acceptedstrings_sam.append(pidStr);
			}
			acceptedstrings_sam.append(".sam");
			acceptedsam.open(acceptedstrings_sam);
			acceptedsam.close();
		}

		if (opts.blastout)
		{
			// blast output
			acceptedstrings_blast = opts.filetype_ar;
			if (opts.pid)
			{
				acceptedstrings_blast.append("_");
				acceptedstrings_blast.append(pidStr);
			}
			acceptedstrings_blast.append(".blast");
			acceptedblast.open(acceptedstrings_blast);
			acceptedblast.close();
		}

		// don't touch the log if only reports are generated
		if (opts.doLog && opts.alirep != Runopts::ALIGN_REPORT::report)
		{
			// statistics file output
			logfile = opts.filetype_ar;
			if (opts.pid)
			{
				logfile.append("_");
				logfile.append(pidStr);
			}
			logfile.append(".log");

			logstream.open(logfile);
			logstream.close();
		}

		if (opts.otumapout)
		{
			// OTU map output file
			ofstream otumap;
			acceptedotumap_file = opts.filetype_ar;
			if (opts.pid)
			{
				acceptedotumap_file.append("_");
				acceptedotumap_file.append(pidStr);
			}
			acceptedotumap_file.append("_otus.txt");
			otumap.open(acceptedotumap_file);
			otumap.close();
		}

		if (opts.de_novo_otu)
		{
			ofstream denovo_otu;
			denovo_otus_file = opts.filetype_ar;
			if (opts.pid)
			{
				denovo_otus_file.append("_");
				denovo_otus_file.append(pidStr);
			}
			denovo_otus_file.append("_denovo.");
			denovo_otus_file.append(readstats.suffix.c_str());

			denovo_otu.open(denovo_otus_file);
			denovo_otu.close();
		}
	}//~if ( ptr_filetype_ar != NULL ) 

	if (opts.filetype_or.size() != 0)
	{
		if (opts.fastxout)
		{
			// output stream for other reads
			ofstream otherreads;
			// add suffix database name to accepted reads file
			if (opts.pid)
			{
				opts.filetype_or += "_";
				opts.filetype_or += pidStr;
			}
			opts.filetype_or += ".";
			opts.filetype_or += readstats.suffix;
			// create the other reads file
			otherreads.open(opts.filetype_or);
			otherreads.close();
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
	const char to_char[5] = { 'A','C','G','T','N' };
	uint32_t id = 0;
	uint32_t mismatches = 0;
	uint32_t gaps = 0;

	// iterate all alignments of the read
	for (int i = 0; i < read.hits_align_info.alignv.size(); ++i)
	{
		uint32_t bitscore = (uint32_t)((float)((refstats.gumbel[refs.num].first)
			* (read.hits_align_info.alignv[i].score1) - log(refstats.gumbel[refs.num].second)) / (float)log(2));

		//double evalue_score = (double)gumbel_K_index_num * full_ref_index_num * full_read_index_num
		// * pow(EXP, (-gumbel_lambda_index_num * result->score1));
		double evalue_score = (double)refstats.gumbel[refs.num].second 
			* refstats.full_ref[refs.num]
			* refstats.full_read[refs.num]
			* std::exp(-refstats.gumbel[refs.num].first * read.hits_align_info.alignv[i].score1);

		std::string refseq = refs.buffer[read.hits_align_info.alignv[i].ref_seq].sequence;
		std::string ref_id = refs.buffer[read.hits_align_info.alignv[i].ref_seq].getId();

		// Blast-like pairwise alignment (only for aligned reads)
		if (opts.blastFormat == BlastFormat::REGULAR) // TODO: global - fix
		{
			acceptedblast << "Sequence ID: ";
			acceptedblast << ref_id; // print only start of the header till first space
			acceptedblast << endl;

			acceptedblast << "Query ID: ";
			acceptedblast << read.getSeqId();
			acceptedblast << endl;

			//fileout << "Score: " << a->score1 << " bits (" << bitscore << ")\t";
			acceptedblast << "Score: " << read.hits_align_info.alignv[i].score1 << " bits (" << bitscore << ")\t";
			acceptedblast.precision(3);
			acceptedblast << "Expect: " << evalue_score << "\t";

			if (read.hits_align_info.alignv[i].strand) acceptedblast << "strand: +\n\n";
			else acceptedblast << "strand: -\n\n";

			if (read.hits_align_info.alignv[i].cigar.size() > 0)
			{
				uint32_t j, c = 0, left = 0, e = 0,
					qb = read.hits_align_info.alignv[i].ref_begin1,
					pb = read.hits_align_info.alignv[i].read_begin1; //mine

				while (e < read.hits_align_info.alignv[i].cigar.size() || left > 0)
				{
					int32_t count = 0;
					int32_t q = qb;
					int32_t p = pb;
					acceptedblast << "Target: ";
					acceptedblast.width(8);
					acceptedblast << q + 1 << "    ";
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
							if (letter == 1) acceptedblast << "-";
							else
							{
								acceptedblast << to_char[(int)refseq[q]];
								++q;
							}
							++count;
							if (count == 60) goto step2;
						}
					}
				step2:
					acceptedblast << "    " << q << "\n";
					acceptedblast.width(20);
					acceptedblast << " ";
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
								if ((char)to_char[(int)refseq[q]] == (char)to_char[(int)read.sequence[p]]) acceptedblast << "|";
								else acceptedblast << "*";
								++q;
								++p;
							}
							else
							{
								acceptedblast << " ";
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
					acceptedblast << "\nQuery: ";
					acceptedblast.width(9);
					acceptedblast << p + 1 << "    ";
					count = 0;
					for (c = e; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
					{
						uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
						uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
						uint32_t l = (count == 0 && left > 0) ? left : length;
						for (j = 0; j < l; ++j)
						{
							if (letter == 2) acceptedblast << "-";
							else
							{
								acceptedblast << (char)to_char[(int)read.sequence[p]];
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
					acceptedblast << "    " << p << "\n\n";
				}
			}
		}
		// Blast tabular m8 + optional columns for CIGAR and query coverage
		else if (opts.blastFormat == BlastFormat::TABULAR)
		{
			// (1) Query
			acceptedblast << read.getSeqId(); // part of the header till first space

			// print null alignment for non-aligned read
			if (opts.print_all_reads && (read.hits_align_info.alignv.size() == 0))
			{
				acceptedblast << "\t*\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0";
				for (uint32_t l = 0; l < opts.blastops.size(); l++)
				{
					if (opts.blastops[l].compare("cigar") == 0)
						acceptedblast << "\t*";
					else if (opts.blastops[l].compare("qcov") == 0)
						acceptedblast << "\t0";
					else if (opts.blastops[l].compare("qstrand") == 0)
						acceptedblast << "\t*";
					acceptedblast << "\n";
				}
				return;
			}

			read.calcMismatchGapId(refs, i, mismatches, gaps, id);
			int32_t total_pos = mismatches + gaps + id;

			acceptedblast << "\t";
			// (2) Subject
			acceptedblast << ref_id << "\t";
			// (3) %id
			acceptedblast.precision(3);
			acceptedblast << (double)id / (mismatches + gaps + id) * 100 << "\t";
			// (4) alignment length
			acceptedblast << (read.hits_align_info.alignv[i].read_end1 - read.hits_align_info.alignv[i].read_begin1 + 1) << "\t";
			// (5) mismatches
			acceptedblast << mismatches << "\t";
			// (6) gap openings
			acceptedblast << gaps << "\t";
			// (7) q.start
			acceptedblast << read.hits_align_info.alignv[i].read_begin1 + 1 << "\t";
			// (8) q.end
			acceptedblast << read.hits_align_info.alignv[i].read_end1 + 1 << "\t";
			// (9) s.start
			acceptedblast << read.hits_align_info.alignv[i].ref_begin1 + 1 << "\t";
			// (10) s.end
			acceptedblast << read.hits_align_info.alignv[i].ref_end1 + 1 << "\t";
			// (11) e-value
			acceptedblast << evalue_score << "\t";
			// (12) bit score
			acceptedblast << bitscore;
			// OPTIONAL columns
			for (uint32_t l = 0; l < opts.blastops.size(); l++)
			{
				// output CIGAR string
				if (opts.blastops[l].compare("cigar") == 0)
				{
					acceptedblast << "\t";
					// masked region at beginning of alignment
					if (read.hits_align_info.alignv[i].read_begin1 != 0) acceptedblast << read.hits_align_info.alignv[i].read_begin1 << "S";
					for (int c = 0; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
					{
						uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
						uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
						acceptedblast << length;
						if (letter == 0) acceptedblast << "M";
						else if (letter == 1) acceptedblast << "I";
						else acceptedblast << "D";
					}

					uint32_t end_mask = read.sequence.length() - read.hits_align_info.alignv[i].read_end1 - 1;
					// output the masked region at end of alignment
					if (end_mask > 0) acceptedblast << end_mask << "S";
				}
				// output % query coverage
				else if (opts.blastops[l].compare("qcov") == 0)
				{
					acceptedblast << "\t";
					acceptedblast.precision(3);
					double coverage = abs(read.hits_align_info.alignv[i].read_end1 - read.hits_align_info.alignv[i].read_begin1 + 1)
						/ read.hits_align_info.alignv[i].readlen;
					acceptedblast << coverage * 100; // (double)align_len / readlen
				}
				// output strand
				else if (opts.blastops[l].compare("qstrand") == 0)
				{
					acceptedblast << "\t";
					if (read.hits_align_info.alignv[i].strand) acceptedblast << "+";
					else acceptedblast << "-";
				}
			}
			acceptedblast << "\n";
		}//~blast tabular m8
	} // ~iterate all alignments
} // ~ Output::report_blast


void Output::writeSamHeader(Runopts & opts)
{
	acceptedsam << "@HD\tVN:1.0\tSO:unsorted\n";

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
	acceptedsam << "@PG\tID:sortmerna\tVN:1.0\tCL:" << opts.cmdline << std::endl;

} // ~Output::writeSamHeader

void Output::report_sam
(
	Runopts & opts,
	References & refs,
	Read & read
)
{
	const char to_char[5] = { 'A','C','G','T','N' };

	// (1) Query
	acceptedsam << read.getSeqId();
	// read did not align, output null string
	if (opts.print_all_reads && (read.hits_align_info.alignv.size() == 0))
	{
		acceptedsam << "\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";
		return;
	}

	// read aligned, output full alignment
	// iterate read alignments
	for (int i = 0; i < read.hits_align_info.alignv.size(); ++i)
	{
		// (2) flag Forward/Reversed
		if (!read.hits_align_info.alignv[i].strand) acceptedsam << "\t16\t";
		else acceptedsam << "\t0\t";
		// (3) Subject
		acceptedsam << refs.buffer[read.hits_align_info.alignv[i].ref_seq].getId();
		// (4) Ref start
		acceptedsam << "\t" << read.hits_align_info.alignv[i].ref_begin1 + 1; // a->ref_begin1
		// (5) mapq
		acceptedsam << "\t" << 255 << "\t";
		// (6) CIGAR
		// output the masked region at beginning of alignment
		if (read.hits_align_info.alignv[i].read_begin1 != 0) 
			acceptedsam << read.hits_align_info.alignv[i].read_begin1 << "S";

		for (int c = 0; c < read.hits_align_info.alignv[i].cigar.size(); ++c)
		{
			uint32_t letter = 0xf & read.hits_align_info.alignv[i].cigar[c];
			uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[i].cigar[c]) >> 4;
			acceptedsam << length;
			if (letter == 0) acceptedsam << "M";
			else if (letter == 1) acceptedsam << "I";
			else acceptedsam << "D";
		}
		uint32_t end_mask = read.sequence.size() - read.hits_align_info.alignv[i].read_end1 - 1; // readlen - a->read_end1
		// output the masked region at end of alignment
		if (end_mask > 0) acceptedsam << end_mask << "S";
		// (7) RNEXT, (8) PNEXT, (9) TLEN
		acceptedsam << "\t*\t0\t0\t";
		// (10) SEQ
		acceptedsam << read.sequence;
		// (11) QUAL
		acceptedsam << "\t";
		// reverse-complement strand
		if (read.quality.size() > 0 && !read.hits_align_info.alignv[i].strand)
		{
			std::reverse(read.quality.begin(), read.quality.end());
			acceptedsam << read.quality;
		}
		else if (read.quality.size() > 0) // forward strand
		{
			acceptedsam << read.quality;
			// FASTA read
		}
		else acceptedsam << "*";

		// (12) OPTIONAL FIELD: SW alignment score generated by aligner
		acceptedsam << "\tAS:i:" << read.hits_align_info.alignv[i].score1;
		// (13) OPTIONAL FIELD: edit distance to the reference
		uint32_t mismatches = 0;
		uint32_t gaps = 0;
		uint32_t id = 0;
		read.calcMismatchGapId(refs, i, mismatches, gaps, id);
		acceptedsam << "\tNM:i:" << mismatches + gaps << "\n";
	} // ~for read.alignments
} // ~Output::report_sam

/* 
 * prototype outputformats.cpp:report_fasta
 *
 * @param reads: 1 or 2 reads (if paired)
 */
void Output::report_fasta(Runopts & opts, std::vector<Read> & reads)
{
	std::stringstream ss;
	double s, f; // for timing

	// output accepted reads
	if ((opts.filetype_ar.size() != 0) && opts.fastxout)
	{
		ss << "    Writing aligned FASTA/FASTQ ... ";
		std::cout << ss.str(); ss.str("");

		TIME(s);
		if (opts.fastxout)
			acceptedreads.open(acceptedstrings, ios::app | ios::binary);

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
					if (acceptedreads.is_open())
					{
						for (Read read: reads)
							acceptedreads << read.header << std::endl << read.sequence << std::endl;
					}
					else
					{
						fprintf(stderr, "  %sERROR%s: [Line %d: %s] file %s could not be opened for writing.\n\n",
							"\033[0;31m", "\033[0m", __LINE__, __FILE__, acceptedstrings.c_str());
						exit(EXIT_FAILURE);
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
					if (acceptedreads.is_open())
					{
						acceptedreads << reads[0].header << std::endl << reads[0].sequence << std::endl;
					}
					else
					{
						fprintf(stderr, "  %sERROR%s: file %s (acceptedstrings) could not be "
							"opened for writing.\n\n", "\033[0;31m", acceptedstrings.c_str(), "\033[0m");
						exit(EXIT_FAILURE);
					}
				}
			} //~if read was accepted
		}//~if not paired-in or paired-out
		if (acceptedreads.is_open()) acceptedreads.close();
		TIME(f);
		ss << " done [" << std::setprecision(2) << (f - s) << " sec]\n"; std::cout << ss.str(); ss.str("");
	}//~if ( ptr_filetype_ar != NULL )

	 // output other reads
	if ((opts.filetype_or.size() != 0) && opts.fastxout)
	{
		ss << "    Writing not-aligned FASTA/FASTQ ... "; std::cout << ss.str(); ss.str("");
		TIME(s);
		otherreads.open(opts.filetype_or, ios::app | ios::binary);
		// pair-ended reads
		if (opts.pairedin || opts.pairedout)
		{
			// neither of the reads were accepted, or exactly one was accepted and pairedout_gv
			if ((!reads[0].hit && !reads[1].hit) ||
				((reads[0].hit ^ reads[1].hit) && opts.pairedout))
			{
				if (otherreads.is_open())
				{
					for (Read read : reads)
						otherreads << read.header << std::endl << read.sequence << std::endl;
				}
				else
				{
					fprintf(stderr, "  %sERROR%s: file %s could not be opened for writing.\n\n", "\033[0;31m", opts.filetype_or.c_str(), "\033[0m");
					exit(EXIT_FAILURE);
				}
			}//~the read was accepted
		}//~if (pairedin_gv || pairedout_gv)
		else // output reads single
		{
			// the read was accepted
			if (!reads[0].hit)
			{
				// accepted reads file output
				if (otherreads.is_open())
				{
					acceptedreads << reads[0].header << std::endl << reads[0].sequence << std::endl;
				}
				else
				{
					fprintf(stderr, "  %sERROR%s: file %s could not be opened for writing.\n\n", "\033[0;31m", opts.filetype_or.c_str(), "\033[0m");
					exit(EXIT_FAILURE);
				}
			}
		} // ~ if (pairedin_gv || pairedout_gv)
		if (otherreads.is_open()) otherreads.close();
		TIME(f);
		ss << " done [" << std::setprecision(2) << (f - s) << " sec]\n"; std::cout << ss.str(); ss.str("");
	}//~if ( opts.fastxout )  
} // ~Output::report_fasta

void Output::report_denovo(Runopts & opts, std::vector<Read> & reads)
{
	std::stringstream ss;
	double s, f; // for timing

	// output reads with < id% alignment (passing E-value) for de novo clustering
	if (denovo_otus_file.size() != 0)
	{
		ss << "    Writing de novo FASTA/FASTQ ... "; std::cout << ss.str(); ss.str("");
		TIME(s);

		denovoreads.open(denovo_otus_file, ios::app | ios::binary);

		// pair-ended reads
		if (opts.pairedin || opts.pairedout)
		{
			// either both reads are accepted, or one is accepted and pairedin_gv
			if (reads[0].hit_denovo || reads[1].hit_denovo && opts.pairedin)
			{
				// output aligned read
				if (denovoreads.is_open())
				{
					for (Read read : reads)
						denovoreads << read.header << std::endl << read.sequence << std::endl;
				}
				else
				{
					fprintf(stderr, "  %sERROR%s: file %s (denovoreads) could not be opened for writing.\n\n", "\033[0;31m", denovo_otus_file.c_str(), "\033[0m");
					exit(EXIT_FAILURE);
				}
			}//~the read was accepted
		}//~if paired-in or paired-out
		else // regular or pair-ended reads don't need to go into the same file
		{
			// the read was accepted
			if (reads[0].hit_denovo)
			{
				// output aligned read
				if (denovoreads.is_open())
				{
					denovoreads << reads[0].header << std::endl << reads[0].sequence << std::endl;
				}
				else
				{
					fprintf(stderr, "  %sERROR%s: file %s (denovoreads) could not be opened for writing.\n\n", "\033[0;31m", denovo_otus_file.c_str(), "\033[0m");
					exit(EXIT_FAILURE);
				}
			} //~if read was accepted
		}//~if not paired-in or paired-out

		if (denovoreads.is_open()) denovoreads.close();

		TIME(f);
		ss << " done [" << std::setprecision(2) << (f - s) << " sec]\n"; std::cout << ss.str(); ss.str("");
	}//~if ( ptr_filetype_ar != NULL )
} // ~Output::report_denovo

void Output::report_biom(){

	biomout.open(biomfile, ios::in);

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
	if (opts.blastout) acceptedblast.open(acceptedstrings_blast);

	if (opts.samout) { 
		acceptedsam.open(acceptedstrings_sam);
		if (!acceptedsam.good())
		{
			fprintf(stderr, "  %sERROR%s: could not open SAM output file for writing.\n", startColor, endColor);
			exit(EXIT_FAILURE);
		}
	}
}

void Output::closefiles()
{
	if (acceptedblast.is_open()) acceptedblast.close();
	if (acceptedsam.is_open()) acceptedsam.close();
}

// called from main. TODO: move into a class?
void generateReports(Runopts & opts)
{
	int N_READ_THREADS = 1;
	int N_PROC_THREADS = 1;
	int loopCount = 0; // counter of total number of processing iterations. TODO: no need here?
	std::stringstream ss;

	ss << "generateReports called. Thread: " << std::this_thread::get_id() << std::endl;
	std::cout << ss.str(); ss.str("");

	ThreadPool tpool(N_READ_THREADS + N_PROC_THREADS);
	KeyValueDatabase kvdb(opts.kvdbPath);
	ReadsQueue readQueue("read_queue", QUEUE_SIZE_MAX, N_READ_THREADS); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", QUEUE_SIZE_MAX, N_PROC_THREADS); // Not used for Reports
	Readstats readstats(opts);
	Output output(opts, readstats);
	Index index;
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
			ss << "Loading index part " << idx_part+1 << "/" << refstats.num_index_parts[index_num] << "  ... ";
			std::cout << ss.str(); ss.str("");
			auto t = std::chrono::high_resolution_clock::now();
			index.load(index_num, idx_part, opts, refstats);
			refs.load(index_num, idx_part, opts, refstats);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t; // ~20 sec Debug/Win
			ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << " sec]\n";
			std::cout << ss.str(); ss.str("");

			for (int i = 0; i < N_READ_THREADS; ++i)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), opts, readQueue, kvdb, loopCount));
			}

			// add processor jobs
			for (int i = 0; i < N_PROC_THREADS; ++i)
			{
				tpool.addJob(ReportProcessor("proc_" + std::to_string(i), readQueue, opts, refs, output, refstats, reportsJob));
			}
			++loopCount;
			tpool.waitAll(); // wait till processing is done on one index part
			index.clear();
			refs.clear();
			writeQueue.reset(N_PROC_THREADS);
		} // ~for(idx_part)
	} // ~for(index_num)
	output.closefiles();
	std::cout << "Done generateReports\n";
} // ~generateReports