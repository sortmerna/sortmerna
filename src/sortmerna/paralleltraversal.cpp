/**
 * @file paralleltraversal.cpp
 * @brief File containing functions for index traversal.
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
 * along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 * @endparblock
 *
 * @contributors Jenya Kopylova   jenya.kopylov@gmail.com
 *               Laurent Noé      laurent.noe@lifl.fr
 *               Pierre Pericard  pierre.pericard@lifl.fr
 *               Daniel McDonald  wasade@gmail.com
 *               Mikaël Salson    mikael.salson@lifl.fr
 *               Hélène Touzet    helene.touzet@lifl.fr
 *               Rob Knight       robknight@ucsd.edu
 */


#include <algorithm>
#include <locale>
#include <iomanip> // output formatting

#include "paralleltraversal.hpp"
#include "load_index.hpp"
#include "kseq.h"
#include "kseq_load.hpp"
#include "traverse_bursttrie.hpp"
#include "alignment.hpp"
#include "alignment2.hpp"
#include "mmap.hpp"

#include "options.hpp"
#include "ThreadPool.hpp"
#include "read.hpp"
#include "readstats.hpp"
#include "index.hpp"
#include "references.hpp"
#include "readsqueue.hpp"
#include "kvdb.hpp"
#include "processor.hpp"
#include "reader.hpp"
#include "writer.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(_WIN32)
#define O_SMR_READ_BIN O_RDONLY | O_BINARY
#else
#define O_SMR_READ_BIN O_RDONLY
#endif

// forward
void writeLog(Runopts & opts, Index & index, Readstats & readstats, Output & output);
void writeAlignments(Runopts & opts, Index & index, References & refs, Readstats & readstats);
void writeAlignmentJob(Index & index, References & refs, Output & output, Readstats & readstats, Read & read); // callback

 // see "heuristic 1" below
 //#define HEURISTIC1_OFF

 /*! @brief Return complement of a nucleotide in
	 integer format.

	 <table>
	  <tr><th>i</th> <th>complement[i]</th></tr>
	  <tr><td>0 (A)</td> <td>3 (T)</td></tr>
	  <tr><td>1 (C)</td> <td>2 (G)</td></tr>
	  <tr><td>2 (G)</td> <td>1 (C)</td></tr>
	  <tr><td>3 (T)</td> <td>0 (A)</td></tr>
	 </table>
  */
char complement[4] = { 3,2,1,0 };

/* @function format_forward()
 * format the forward read into a string on same alphabet without '\n'
 *
 * */
void format_forward(char* read_seq, char* myread, char filesig)
{
	// FASTA
	if (filesig == '>')
	{
		while ((*read_seq != '\0') && (*read_seq != '>'))
		{
			if (*read_seq != '\n' || *read_seq != '\r') *myread++ = nt_table[(int)*read_seq];
			read_seq++;
		}
		*myread = '\n'; // end of read marked by newline
	}
	// FASTQ
	else
	{
		while (*read_seq != '\n' &&  *read_seq != '\r') { *myread++ = nt_table[(int)*read_seq++]; }
		*myread = '\n'; //end of read marked by newline
	}
}

/* @function format_rev()
 * format the reverse-complement read into a string without '\n'
 *
 * */
void format_rev(char* start_read, char* end_read, char* myread, char filesig)
{
	int8_t rc_table[128] = {
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 3, 4, 2, 4, 4, 4, 1, 4, 4, 4, 4, 4, 4, 4, 4,
	  4, 4, 4, 4, 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
	};

	// FASTA
	if (filesig == '>')
	{
		while (end_read != start_read)
		{
			if (*end_read != '\n' && *end_read != '\r') *myread++ = rc_table[(int)*end_read];
			end_read--;
		}
		*myread++ = rc_table[(int)*end_read];
		*myread = '\n';
	}
	// FASTQ
	else
	{
		while (*end_read != '\n' && *end_read != '\r') { *myread++ = rc_table[(int)*end_read--]; }
		*myread = '\n';
	}
}

/* TODO: remove -> Readstats::check_file_format
 *
 * @function check_file_format()
 * see documentation in paralleltraversal.hpp
 *
 * */
bool check_file_format(const char* inputreads, char& filesig)
{
	bool exit_early = false;
#ifdef HAVE_LIBZ
	// Check file format (if ZLIB supported)
	gzFile fp = gzopen(inputreads, "r");
	kseq_t *seq = kseq_init(fp);
#else
	FILE* fp = fopen(inputreads, "rb");
	kseq_t *seq = kseq_init(fileno(fp));
#endif
	int l;
	if ((l = kseq_read(seq)) >= 0)
		filesig = seq->last_char;
	else
	{
		fprintf(stderr, "  %sERROR%s: Line %d: %s unrecognized file format or empty file %s\n\n",
			startColor, "\033[0m", __LINE__, __FILE__, inputreads);
		exit_early = true;
	}
	kseq_destroy(seq);
#ifdef HAVE_LIBZ
	gzclose(fp);
#else
	fclose(fp);
#endif
	return exit_early;
}//~check_file_format()


/* @function compute_read_stats()
 * see documentation in paralleltraversal.hpp
 *
 * */
void compute_read_stats(
	char* inputreads,
	uint64_t& number_total_read,
	uint64_t& full_read_main,
	off_t& full_file_size
) {
#ifdef HAVE_LIBZ
	// Count total number of reads and their combined length
	// (if ZLIB is supported)
	gzFile fp = gzopen(inputreads, "r");
	kseq_t *seq = kseq_init(fp);
#else
	// Count total number of reads and their combined length
	// (if ZLIB is not supported)
	FILE* fp = fopen(inputreads, "rb");
	kseq_t *seq = kseq_init(fileno(fp));
#endif
	int l;
	while ((l = kseq_read(seq)) >= 0) {
		full_read_main += seq->seq.l; // seq.l is sequence length
		number_total_read++;
		// compute size of all reads to store in memory
		// + 7 (4 possible new lines, 2 symbols > or @ and +, space for comment)
		if (!map_size_set_gv)
			full_file_size += (seq->name.l + seq->comment.l + seq->seq.l + seq->qual.l + 7);
	}
	if (l == -2)
	{
		fprintf(stderr, "  %sERROR%s: Line %d: %s could not read reads file - %s\n\n",
			startColor, "\033[0m", __LINE__, __FILE__, strerror(errno));
		exit(EXIT_FAILURE);
	}
	kseq_destroy(seq);
#ifdef HAVE_LIBZ
	gzclose(fp);
#else
	fclose(fp);
#endif
}//~compute_read_stats()


// Callback run in a Processor thread
void parallelTraversalJob(Index & index, References & refs, Output & output, Readstats & readstats, Read & read)
{
	// for reverse reads
	if (!forward_gv)
	{
		// output the first num_alignments_gv alignments
		if (num_alignments_gv > 0)
		{
			// all num_alignments_gv alignments have been output
			if (read.num_alignments < 0) return;
		}
		// the maximum scoring alignment has been found, go to next read
		// (unless all alignments are being output)
		else if ((num_best_hits_gv > 0) &&	(min_lis_gv > 0) 
			&& (read.max_SW_score == num_best_hits_gv))
			return;
	}

// TODO: Delete --------------------------------------------------->
	//char myread[READLEN] = ""; // sequence converted to integer alphabet (0,1,2,3)  read.isequence
	//const char* str = read.sequence.c_str();
	//char* _myread = &myread[0];
	// keep record of what position on the read an ambiguous char exists,
	// the char at this position will be changed to '4' during
	// sequence alignment
	//int32_t ambiguous_nt[READLEN] = { 0 }; // Read::ambiguous_nt - vector
	// size of ambiguous_nt
	//uint32_t size_ambiguous_nt = 0;
	// flag to count only 1 alignment per read
	bool read_to_count = true; // passed directly to compute_lis_alignment. What's the point?
	// length of read
	//uint32_t readlen = 0;

	// change the read into an integer alphabet -- FASTA
	{
		// block deleted. Moved to Reads::isequence + Reads::seqToIntStr
	}
// <---------------------------------------------------------------- Delete

	// find the minimum sequence length
	if (read.sequence.size() < readstats.min_read_len) readstats.min_read_len = static_cast<uint32_t>(read.sequence.size());
	// find the maximum sequence length
	if (read.sequence.size()  > readstats.max_read_len) readstats.max_read_len = static_cast<uint32_t>(read.sequence.size());

	// the read length is too short
	if (read.sequence.size()  < index.lnwin[index.index_num])
	{
		std::stringstream ss;
		ss << std::this_thread::get_id();
		fprintf(stderr, "\n  %sWARNING%s: Processor thread: %s The read %d is shorter "
			"than %u nucleotides, by default it will not be searched\n", "\033[0;33m", "\033[0m", 
			ss.str().c_str(), read.id, index.lnwin[index.index_num]);
		read.isValid = false;
		return;
	}

	// create the reverse complement of the sequence
	if (!forward_gv)
	{
		if (!read.reversed) read.revIntStr();
	}//~if (REVERSE)

	 // array of positions of window hits on the reference sequence
	//vector< id_win > id_win_hits; // Read::id_win_hits
	// number of windows hit to the reference sequence(s)
	//uint32_t readhit = 0; read.readhit
	uint32_t windowshift = index.opts.skiplengths[index.index_num][0];
	// keep track of windows which have been already traversed in the burst trie
	vector<bool> read_index_hits(read.sequence.size());
	
	uint32_t pass_n = 0; // Pass number (possible value 0,1,2)
	// the maximum SW score attainable for this read
	uint32_t max_SW_score = read.sequence.size() *index.opts.match;

	std::vector<MYBITSET> vbitwindowsf;
	std::vector<MYBITSET> vbitwindowsr;

	// TODO: these 2 lines needs to be calculated once per index part. Move to index or some new calculation object
	uint32_t bit_vector_size = (index.partialwin[index.index_num] - 2) << 2; // on each index part. Init in Main, accessed in Worker
	uint32_t offset = (index.partialwin[index.index_num] - 3) << 2; // on each index part

	uint32_t minoccur = 0; // TODO: never updated. Always 0. What's the point?

	// loop for each new Pass to granulate seed search intervals
	for	(bool search = true; search; )
	{  
		uint32_t numwin = (read.sequence.size() 
			- index.lnwin[index.index_num] 
			+ windowshift) / windowshift; // number of k-mer windows fit along the sequence

		uint32_t win_index = 0; // index of the window's first char in the sequence e.g. 0, 18, 36 if window.length = 18
		// iterate over windows of the template string
		for (uint32_t win_num = 0; win_num < numwin; win_num++)
		{
			// skip position, seed at this position has already been searched for in a previous Pass
			if (read_index_hits[win_index]) goto check_score;
			// search position, set search bit to true
			else read_index_hits[win_index].flip();
			{
				// this flag it set to true if a match is found during
				// subsearch 1(a), to skip subsearch 1(b)
				bool accept_zero_kmer = false;
				// ids for k-mers that hit the database
				vector< id_win > id_hits; // TODO: why not to add directly to 'id_win_hits'?
				vbitwindowsf.resize(bit_vector_size);
				std::fill(vbitwindowsf.begin(), vbitwindowsf.end(), 0);

				init_win_f(&read.isequence[win_index + index.partialwin[index.index_num]],
					&vbitwindowsf[0],
					&vbitwindowsf[4],
					index.numbvs[index.index_num]);

				uint32_t keyf = 0;
				char *keyf_ptr = &read.isequence[win_index];
				// build hash for first half windows (foward and reverse)
				// hash is just a numeric value formed by the chars of a string consisting of '0','1','2','3'
				// e.g. "2233012" -> b10.1011.1100.0110 = x2BC6 = 11206
				for (uint32_t g = 0; g < index.partialwin[index.index_num]; g++)
				{
					(keyf <<= 2) |= (uint32_t)(*keyf_ptr);
					++keyf_ptr;
					//(keyf <<= 2) |= (uint32_t)*keyf_ptr++; // TODO: How did this work? And it did!
				}
				// do traversal if the exact half window exists in the burst trie
				if ((index.lookup_tbl[keyf].count > minoccur) && (index.lookup_tbl[keyf].trie_F != NULL))
				{
					/* subsearch (1)(a) d([p_1],[w_1]) = 0 and d([p_2],[w_2]) <= 1;
					*
					*  w = |------ [w_1] ------|------ [w_2] ------|
					*  p = |------ [p_1] ------|------ [p_2] ----| (0/1 deletion in [p_2])
					*              or
					*    = |------ [p_1] ------|------ [p_2] ------| (0/1 match/substitution in [p_2])
					*        or
					*    = |------ [p_1] ------|------ [p_2] --------| (0/1 insertion in [p_2])
					*
					*/   
					traversetrie_align(
						index.lookup_tbl[keyf].trie_F,
						0,
						0,
						// win2f_k1_ptr
						&vbitwindowsf[0],
						// win2f_k1_full
						&vbitwindowsf[offset],
						accept_zero_kmer,
						id_hits,
						read.id,
						win_index,
						index.partialwin[index.index_num]
					);        
				}//~if exact half window exists in the burst trie

				 // only search if an exact match has not been found
				if (!accept_zero_kmer)
				{
					vbitwindowsr.resize(bit_vector_size);
					std::fill(vbitwindowsr.begin(), vbitwindowsr.end(), 0);
					// build the first bitvector window
					init_win_r(&read.isequence[win_index + index.partialwin[index.index_num] - 1],
						&vbitwindowsr[0],
						&vbitwindowsr[4],
						index.numbvs[index.index_num]);
					uint32_t keyr = 0;
					char *keyr_ptr = &read.isequence[win_index + index.partialwin[index.index_num]];
					// build hash for first half windows (foward and reverse)
					for (uint32_t g = 0; g < index.partialwin[index.index_num]; g++)
					{
						(keyr <<= 2) |= (uint32_t)(*keyr_ptr-'0');
						++keyr_ptr;
						//(keyr <<= 2) |= (uint32_t)*keyr_ptr++;
					}
					// continue subsearch (1)(b)
					if ((index.lookup_tbl[keyr].count > minoccur) && (index.lookup_tbl[keyr].trie_R != NULL))
					{
						/* subsearch (1)(b) d([p_1],[w_1]) = 1 and d([p_2],[w_2]) = 0;
						*
						*  w =    |------ [w_1] ------|------ [w_2] -------|
						*  p =      |------- [p_1] ---|--------- [p_2] ----| (1 deletion in [p_1])
						*              or
						*    =    |------ [p_1] ------|------ [p_2] -------| (1 match/substitution in [p_1])
						*        or
						*    = |------- [p_1] --------|---- [p_2] ---------| (1 insertion in [p_1])
						*
						*/
						traversetrie_align(
							index.lookup_tbl[keyr].trie_R,
							0,
							0,
							/* win1r_k1_ptr */
							&vbitwindowsr[0],
							/* win1r_k1_full */
							&vbitwindowsr[offset],
							accept_zero_kmer,
							id_hits,
							read.id,
							win_index,
							index.partialwin[index.index_num]);
					}//~if exact half window exists in the reverse burst trie                    
				}//~if (!accept_zero_kmer)                         
				 // associate the ids with the read window number
				if (!id_hits.empty())
				{
					for (uint32_t i = 0; i < id_hits.size(); i++)
					{
						read.id_win_hits.push_back(id_hits[i]);
					}
					read.readhit++;
				}
			}
		check_score:
			// continue read analysis if threshold seeds were matched
			if (win_num == numwin - 1)
			{
				compute_lis_alignment2(
					read, index, refs, readstats, output,
					//size_ambiguous_nt, TODO: remove. Not used no more.
					// readhit, Read::readhit
					//id_win_hits,  Read::id_win_hits
					//positions_tbl, Index
					//read_max_SW_score, Read::max_SW_score
					search, // returns False if the alignment is found -> stop searching
					//best_x, Read::best
					// readn,  Read::id
					//num_alignments_x,  Read::num_alignments
					//readlen,  Read::sequence.length
					//lnwin[index_num],  Index::lnwin
					//index_num,  Index::index_num
					//reference_seq_len, // References
					//myread,  Read::isequence
					//ambiguous_nt,   Read::ambiguous_nt vector
					//scoring_matrix, Read::scoring_matrix
					//reference_seq,  References
					//gap_open, Runopts::gap_open
					//gap_extension, Runops
					//minimal_score[index_num], Index
					//read_hits,  Read::hit
					//total_reads_mapped,  Output::
					//reads_matched_per_db, Output::
					//part, Index::part
					//read_hits_align_info, Read::hits_align_info
					max_SW_score,
					read_to_count
					//total_reads_mapped_cov, Output::
					//read_hits_denovo,  Read::hit_denovo
					//filesig, Read::
					//strs,  total number of reads - no more
					//file_s,  Mmap related
					//file_sections, Mmap related
					//reads, Read
					//finalnt,  marks the end of mmap'd region
					//gumbel[index_num].first,  Index::gambel[index_num].first
					//gumbel[index_num].second, Index::gambel[index_num].second
					//full_ref[index_num], Index::
					//full_read[index_num], Index::
					//acceptedblast,  Output::acceptedblast
					//acceptedsam  OUtput::acceptedsam
				);
				// the read was not accepted at current window skip length,
				// decrease the window skip length
				if (search)
				{
					// last (3rd) Pass has been made
					if (pass_n == 2) search = false;
					else
					{
						// the next interval size equals to the current one, skip it
						while ((pass_n < 3) &&
							(index.opts.skiplengths[index.index_num][pass_n] == index.opts.skiplengths[index.index_num][pass_n + 1])) ++pass_n;
						if (++pass_n > 2) search = false;
						// set interval skip length for next Pass
						else windowshift = index.opts.skiplengths[index.index_num][pass_n];
					}
				}
				break; // do not offset final window on read
			}//~( win_num == NUMWIN-1 )
			win_index += windowshift;
		}//~for (each window)                
		 //~while all three window skip lengths have not been tested, or a match has not been found
	}// while (search);

	// the read didn't align (for --num_alignments [INT] option),
	// output null alignment string
	if (!read.hit && !forward_gv && (num_alignments_gv > -1))
	{
		// do not output read for de novo OTU clustering
		// (it did not pass the E-value threshold)
		if (de_novo_otu_gv && read.hit_denovo) read.hit_denovo = !read.hit_denovo; // flip
		// output null alignment string
		if (print_all_reads_gv)
		{
#pragma omp critical
			{
				s_align* null_alignment = NULL;
				if (index.opts.blastout && blast_tabular)
				{
					report_blast(output.acceptedblast, // blast output file
						null_alignment, // SW alignment cigar
						read.header.c_str(), //reads[readn - 1] + 1, //read name
						0, // read sequence (in integer format)
						0, // read quality
						0, // reference name
						0, // reference sequence
						0, // e-value score
						0, // read length (to compute the masked regions)
						0, // bitscore
						0, // forward or reverse strand
						0, // %id
						0, // %query coverage
						0, // number of mismatches
						0); // number of gaps
				}
				if (index.opts.samout)
				{
					report_sam(output.acceptedsam, // sam output file
						null_alignment, // SW alignment cigar
						read.header.c_str(), //reads[readn - 1] + 1, // read name
						0, // read sequence (in integer format)
						0, // read quality
						0, // reference name
						0, // reference sequence
						0, // read length (to compute the masked regions)
						0, // forward or reverse strand
						0); // edit distance
				}
			}//~if print_all_reads_gv
		}// allow writing to file 1 thread at a time
	}//~if read didn't align

	if (de_novo_otu_gv && read.hit_denovo)
		++readstats.total_reads_denovo_clustering;
} // ~parallelTraversalJob


void paralleltraversal2(Runopts & opts)
{
	unsigned int numCores = std::thread::hardware_concurrency(); // find number of CPU cores
	std::cout << "CPU cores on this machine: " << numCores << std::endl; // 8

	// Init thread pool with the given number of threads
	int numThreads = 2 * opts.num_fread_threads + opts.num_proc_threads;
	if (numThreads > numCores)
		printf("WARN: Number of cores: %d is less than number allocated threads %d", numCores, numThreads);

	ThreadPool tpool(numThreads);
	KeyValueDatabase kvdb(opts.kvdbPath);
	ReadsQueue readQueue("read_queue", QUEUE_SIZE_MAX, 1); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", QUEUE_SIZE_MAX, opts.num_proc_threads); // shared: Processor pushes, Writer pops
	Readstats readstats(opts);
	Output output(opts, readstats);
	Index index(opts, readstats, output);
	References refs(opts, index);

	int loopCount = 0; // counter of total number of processing iterations

	// perform alignment
	// loop through every index passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); index_num++)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < index.num_index_parts[index_num]; idx_part++)
		{
			eprintf("\tLoading index part %d/%u ... ", idx_part + 1, index.num_index_parts[index_num]);
			auto t = std::chrono::high_resolution_clock::now();
			index.load(index_num, idx_part);
			refs.load(index_num, idx_part);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t; // ~20 sec Debug/Win
//			std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - t);
			eprintf("done [%.2f sec]\n", elapsed.count());

			for (int i = 0; i < opts.num_fread_threads; i++)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), opts, readQueue, opts.readsfile, kvdb, loopCount));
				tpool.addJob(Writer("writer_" + std::to_string(i), writeQueue, kvdb));
			}

			// add processor jobs
			for (int i = 0; i < opts.num_proc_threads; i++)
			{
				tpool.addJob(Processor("proc_" + std::to_string(i), readQueue, writeQueue, readstats, index, refs, output, parallelTraversalJob));
			}
			++loopCount;
		} // ~for(idx_part)
	} // ~for(index_num)

	tpool.waitAll(); // wait till processing is done

	writeLog(opts, index, readstats, output);

#if 0
	// Write reports
	// loop through every index passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); index_num++)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < index.num_index_parts[index_num]; idx_part++)
		{
			eprintf("\tLoading index part %d/%u ... ", idx_part + 1, index.num_index_parts[index_num]);
			auto t = std::chrono::high_resolution_clock::now();
			index.load(index_num, idx_part);
			refs.load(index_num, idx_part);
			std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - t; // ~20 sec Debug/Win
			eprintf("done [%.2f sec]\n", elapsed.count());

			for (int i = 0; i < opts.num_fread_threads; i++)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), readQueue, opts.readsfile, kvdb, loopCount));
			}

			// add processor jobs
			for (int i = 0; i < opts.num_proc_threads; i++)
			{
				tpool.addJob(Processor("proc_" + std::to_string(i), readQueue, writeQueue, readstats, index, refs, output, writeAlignmentJob));
			}
			++loopCount;
		} // ~for(idx_part)
	} // ~for(index_num)
#endif
} // ~paralleltraversal2

/*! @fn paralleltraversal() */
void paralleltraversal(
	char* inputreads,
	bool have_reads_gz,
	char* ptr_filetype_ar, // aligned reads file
	char* ptr_filetype_or, // rejected reads file
	long match,
	long mismatch,
	long gap_open,
	long gap_extension,
	long score_N,
	vector<vector<uint32_t>>& skiplengths,
	int argc,
	char **argv,
	bool yes_SQ,
	vector<pair<string, string>>& myfiles,
	bool exit_early)
{
	// the offset from the start of the reads file for mmap
	off_t offset_map = 0;
	// the size of the full reads file (in bytes)
	off_t full_file_size = 0;
	// total number of nucleotides in all reads
	uint64_t full_read_main = 0;
	// total number of reads
	uint64_t number_total_read = 0;
	// total number of reads mapped passing E-value threshold
	uint64_t total_reads_mapped = 0;
	// total number of reads mapped passing E-value threshold and
	// %id and/or %query coverage thresholds
	uint64_t total_reads_mapped_cov = 0;
	// total number of reads for de novo clustering
	uint64_t total_reads_denovo_clustering = 0;
	// the minimum occurrences of a (L/2)-mer required to allow
	// search for a seed of length L in the burst tries
	uint32_t minoccur = 0; // TODO: never updated i.e. always 0. What's the point?
	// for timing different processes
	double s, f;
	// the comparing character used in parsing the reads file
	char filesig;
	// minimum read length for log statistics
	uint32_t min_read_len = READLEN;
	// maximum read length for log statistics
	uint32_t max_read_len = 0;
	// mean read length for log statistics
	uint32_t mean_read_len = 0;
	// the size of the sliding window on the full file, ~1GB
	off_t partial_file_size = 0;
	// size of the remainder of the file (last window) which
	// is less than pow(2,30) bytes
	off_t last_part_size = 0;
	// number of file sections to mmap
	uint32_t file_sections = 0;
	// index for file_sections
	uint32_t file_s = 0;
	exit_early = check_file_format(inputreads, filesig);

	// File format supported (FASTA or FASTQ), continue
	if (!exit_early)
	{
		eprintf("\n  Computing read file statistics ...");
		TIME(s);
		compute_read_stats(inputreads,
			number_total_read,
			full_read_main,
			full_file_size);

		// find the mean sequence length
		mean_read_len = full_read_main / number_total_read;

		// check there are an even number of reads for --paired-in
		// and --paired-out options to work
		if ((number_total_read % 2 != 0) && (pairedin_gv || pairedout_gv))
		{
			fprintf(stderr, "\n    %sWARNING%s: for --paired-in and --paired-out options, the number of reads must be even.\n",
				"\033[0;33m", "\033[0m");
			fprintf(stderr, "    There are %lld reads in your file.\n", number_total_read);
			fprintf(stderr, "    Reads will still be processed and output, but paired-reads may be split.\n\n");
			pairedin_gv = false;
			pairedout_gv = false;
		}

		// setup for mmap
		if (map_size_set_gv)
		{
			// full_file_size for mmap is the exact file size
			int fd = open(inputreads, O_SMR_READ_BIN);
			// find the size of the total file
			if ((full_file_size = lseek(fd, 0L, SEEK_END)) == -1)
			{
				fprintf(stderr, "  %sERROR%s: Line %d: %s Could not seek reads file!\n\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			if (lseek(fd, 0L, SEEK_SET) == -1)
			{
				fprintf(stderr, "  %sERROR%s: Line %d: %s Could not seek set the reads file!\n\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			close(fd);
			partial_file_size = full_file_size;
			last_part_size = full_file_size%map_size_gv;
			// if the full_file_size is bigger than m*PAGE_SIZE, mmap
			// the file by 'windows' of size partial_file_size,
			// otherwise keep the full_file_size
			if ((file_sections = ceil((double)full_file_size / (double)(map_size_gv))) > 1)
				partial_file_size = map_size_gv;
		}
		else
			file_sections = 1;
		TIME(f);
		eprintf("  done [%.2f sec]\n", (f - s));
		if (map_size_set_gv)
		{
			eprintf("  size of reads file: %lu bytes\n", (unsigned long int)full_file_size);
			eprintf("  partial section(s) to be executed: %d of size %lu bytes \n",
				file_sections, (unsigned long int)partial_file_size);
		}
	}//~if (!exit_early)

// Output -------------------------------------------------------------------->
	// output streams for aligned reads (FASTA/FASTQ, SAM and BLAST-like)
	ofstream acceptedreads; // TODO: move to Output
	ofstream acceptedsam; // TODO: move to Output
	ofstream acceptedblast; // TODO: move to Output

	// determine the suffix (fasta, fastq, ...) of aligned strings
	char suffix[20] = "out";
	char *suff = strrchr(inputreads, '.');
	if (suff != NULL && !have_reads_gz)
		strcpy(suffix, suff + 1);
	else if (filesig == '>')
		strcpy(suffix, "fasta");
	else
		strcpy(suffix, "fastq");
	suff = NULL;

	char *acceptedstrings = NULL; // --> Output
	char *acceptedstrings_sam = NULL; // --> Output
	char *acceptedstrings_blast = NULL; // --> Output
	char *logoutfile = NULL; // --> Output
	char *denovo_otus_file = NULL; // --> Output
	char *acceptedotumap_file = NULL; // --> Output

	// attach pid to output files
	char pidStr[4000];
	if (pid_gv)
	{
		int32_t pid = getpid();
		sprintf(pidStr, "%d", pid);
	}

	// associate the streams with reference sequence file names
	if (ptr_filetype_ar != NULL)
	{
		if (fastxout_gv)
		{
			// fasta/fastq output
			acceptedstrings = new char[1000]();
			if (acceptedstrings == NULL)
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for acceptedstrings\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			strcpy(acceptedstrings, ptr_filetype_ar);
			if (pid_gv)
			{
				strcat(acceptedstrings, "_");
				strcat(acceptedstrings, pidStr);
			}
			strcat(acceptedstrings, ".");
			strcat(acceptedstrings, suffix);

			acceptedreads.open(acceptedstrings);
			acceptedreads.close();
		}

		if (samout_gv)
		{
			// sam output
			acceptedstrings_sam = new char[1000]();
			if (acceptedstrings_sam == NULL)
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for acceptedstrings_sam\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			strcpy(acceptedstrings_sam, ptr_filetype_ar);
			if (pid_gv)
			{
				strcat(acceptedstrings_sam, "_");
				strcat(acceptedstrings_sam, pidStr);
			}
			strcat(acceptedstrings_sam, ".sam");
			acceptedsam.open(acceptedstrings_sam);
			acceptedsam.close();
		}

		if (blastout_gv)
		{
			// blast output
			acceptedstrings_blast = new char[1000]();
			if (acceptedstrings_blast == NULL)
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for acceptedstrings_blast\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			strcpy(acceptedstrings_blast, ptr_filetype_ar);
			if (pid_gv)
			{
				strcat(acceptedstrings_blast, "_");
				strcat(acceptedstrings_blast, pidStr);
			}
			strcat(acceptedstrings_blast, ".blast");
			acceptedblast.open(acceptedstrings_blast);
			acceptedblast.close();
		}

		if (logout_gv)
		{
			// statistics file output
			ofstream logstream;
			logoutfile = new char[1000]();
			if (logoutfile == NULL)
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for acceptedstrings_blast\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			strcpy(logoutfile, ptr_filetype_ar);
			if (pid_gv)
			{
				strcat(logoutfile, "_");
				strcat(logoutfile, pidStr);
			}
			strcat(logoutfile, ".log");

			logstream.open(logoutfile);
			logstream.close();
		}

		if (otumapout_gv)
		{
			// OTU map output file
			ofstream otumap;
			acceptedotumap_file = new char[1000]();
			if (acceptedotumap_file == NULL)
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for acceptedotumap\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			strcpy(acceptedotumap_file, ptr_filetype_ar);
			if (pid_gv)
			{
				strcat(acceptedotumap_file, "_");
				strcat(acceptedotumap_file, pidStr);
			}
			strcat(acceptedotumap_file, "_otus.txt");
			otumap.open(acceptedotumap_file);
			otumap.close();
		}

		if (de_novo_otu_gv)
		{
			ofstream denovo_otu;
			denovo_otus_file = new char[1000]();
			if (denovo_otus_file == NULL)
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for denovo_otus_file_name\n",
					startColor, "\033[0m", __LINE__, __FILE__);
				exit(EXIT_FAILURE);
			}
			strcpy(denovo_otus_file, ptr_filetype_ar);
			if (pid_gv)
			{
				strcat(denovo_otus_file, "_");
				strcat(denovo_otus_file, pidStr);
			}
			strcat(denovo_otus_file, "_denovo.");
			strcat(denovo_otus_file, suffix);

			denovo_otu.open(denovo_otus_file);
			denovo_otu.close();
		}
	}//~if ( ptr_filetype_ar != NULL ) 

	if (ptr_filetype_or != NULL)
	{
		if (fastxout_gv)
		{
			// output stream for other reads
			ofstream otherreads;
			// add suffix database name to accepted reads file
			if (pid_gv)
			{
				strcat(ptr_filetype_or, "_");
				strcat(ptr_filetype_or, pidStr);
			}
			strcat(ptr_filetype_or, ".");
			strcat(ptr_filetype_or, suffix);
			// create the other reads file
			otherreads.open(ptr_filetype_or);
			otherreads.close();
		}
	}
// <----------------------------------------------------------- Output

	// empty output files created, exit program
	if (exit_early)
	{
		fprintf(stderr, "  The input reads file or reference file is empty, "
			"or the reads file is not in FASTA or FASTQ format, "
			"no analysis could be made.\n");
		// output parameters to log file
		if ((ptr_filetype_ar != NULL) && logout_gv)
		{
			FILE* bilan = fopen(logoutfile, "w");
			if (bilan == NULL)
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not open file %s \n", startColor, "\033[0m",
					__LINE__, __FILE__, logoutfile);
				exit(EXIT_FAILURE);
			}
			fprintf(bilan, "  The input reads file or reference file is empty, "
				"or the reads file is not in FASTA or FASTQ format, "
				"no analysis could be made.\n");
			fclose(bilan);
		}
		exit(EXIT_SUCCESS);
	}

	// Move to Read::scoring_matrix
	// created/destroied in main thread. Accessed in Worker threads.
	int8_t* scoring_matrix = (int8_t*)calloc(25, sizeof(int8_t));
	{
		int32_t l, k, m;
		// initialize Smith-Waterman scoring matrix for genome sequences
		for (l = k = 0; l < 4; ++l)
		{
			for (m = 0; m < 4; ++m)
				scoring_matrix[k++] = l == m ? match : mismatch; // weight_match : weight_mismatch (must be negative)
			scoring_matrix[k++] = score_N; // ambiguous base
		}
		for (m = 0; m < 5; ++m) scoring_matrix[k++] = score_N; // ambiguous base
	}

	// the number of parts an index was divided into to fit into specified memory,
	// for each reference database searched
	uint32_t num_databases = myfiles.size(); // Index
	// number of index parts per database
	vector<uint16_t> num_index_parts(num_databases, 0); // Index
	// length correction for reference database per database
	vector<uint64_t> full_ref(num_databases, 0); // Index
	// length correction for all reads per database
	vector<uint64_t> full_read(num_databases, full_read_main); // Index
	// L-mer length per database
	vector<uint32_t> lnwin(num_databases, 0); // Index
	// L/2-mer length per database
	vector<uint32_t> partialwin(num_databases, 0); // Index
	// minimal SW score in order to reach threshold E-value
	vector<uint32_t> minimal_score(num_databases, 0); // Index
	// array of structs storing information on which sequences from the original FASTA file were added to each index part
	vector<vector<index_parts_stats> > index_parts_stats_vec; // Index
	// Gumbel parameters lambda and K, respectively
	vector<pair<double, double> > gumbel(num_databases, pair<double, double>(-1.0, -1.0)); // Index
	// total number of full bitvectors in [w_1] reverse and [w_2] forward
	vector<uint64_t> numbvs(num_databases, 0); // Index
	// total number of reads matched for each database in list --ref
	vector<uint64_t> reads_matched_per_db(num_databases, 0); // Index
	// number of reference sequences in each index FASTA file
	vector<uint64_t> numseq(num_databases, 0); // Index
	// set the same skiplengths for all reference files (if the user uses option --passes)
	if (skiplengths.empty())
	{
		vector<uint32_t> skiplengths_v(3, 0);
		skiplengths.push_back(skiplengths_v);
	}
	for (uint32_t i = 0; i < myfiles.size() - 1; i++)
		skiplengths.push_back(skiplengths[0]);

	// TODO: moved to Index::load_stats
	// add header lines to SAM output file.
	load_index_stats(myfiles,
		argv,
		argc,
		yes_SQ,
		acceptedstrings_sam,
		(long)match,
		(long)mismatch,
		(long)gap_open,
		(long)gap_extension,
		skiplengths,
		num_index_parts,
		index_parts_stats_vec,
		full_ref,
		full_read,
		lnwin,
		partialwin,
		minimal_score,
		number_total_read,
		gumbel,
		numbvs,
		numseq);

	// some info on chosen parameters
	eprintf("  Parameters summary:\n");
	eprintf("    Number of seeds = %d\n", seed_hits_gv);
	eprintf("    Edges = %d", edges_gv);
	if (as_percent_gv)
		eprintf(" (as percent)\n");
	else
		eprintf(" (as integer)\n");
	eprintf("    SW match = %ld\n", match);
	eprintf("    SW mismatch = %ld\n", mismatch);
	eprintf("    SW gap open penalty = %ld\n", gap_open);
	eprintf("    SW gap extend penalty = %ld\n", gap_extension);
	eprintf("    SW ambiguous nucleotide = %ld", score_N);
	if (score_N > 0) eprintf(" %sWarning!%s Positive score set for ambiguous nucleotides.\n", "\033[0;33m", "\033[0m");
	else eprintf("\n");
	if (yes_SQ)
		eprintf("    SQ tags are output\n");
	else
		eprintf("    SQ tags are not output\n");
#ifdef _OPENMP
	eprintf("    Number of threads = %d\n", numcpu_gv);
#else
	eprintf("    Number of threads = 1 (OpenMP is not supported with your current C++ compiler).\n");
#endif
	if (pid_gv)
	{
		eprintf("    Current process pid = %d\n", getpid());
	}

	// output parameters to log file
	if ((ptr_filetype_ar != NULL) && logout_gv)
	{
		FILE* bilan = fopen(logoutfile, "w");
		if (bilan == NULL)
		{
			fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not open file %s \n",
				startColor, "\033[0m", __LINE__, __FILE__, logoutfile);
			exit(EXIT_FAILURE);
		}
		time_t q = time(0);
		struct tm * now = localtime(&q);
		fprintf(bilan, " %s\n", asctime(now));
		fprintf(bilan, " Command: ");
		for (int j = 0; j < argc; j++) fprintf(bilan, "%s ", argv[j]);
		fprintf(bilan, "\n");
		// some info on chosen parameters
		fprintf(bilan, " Process pid = %d\n", getpid());
		fprintf(bilan, " Parameters summary:\n");
		for (uint32_t index_num = 0; index_num < myfiles.size(); index_num++)
		{
			fprintf(bilan, "    Index: %s\n", (char*)(myfiles[index_num].second).c_str());
			fprintf(bilan, "     Seed length = %d\n", lnwin[index_num]);
			fprintf(bilan, "     Pass 1 = %d, Pass 2 = %d, Pass 3 = %d\n", skiplengths[index_num][0], skiplengths[index_num][1], skiplengths[index_num][2]);
			fprintf(bilan, "     Gumbel lambda = %f\n", (gumbel[index_num].first));
			fprintf(bilan, "     Gumbel K = %f\n", (gumbel[index_num].second));
			fprintf(bilan, "     Minimal SW score based on E-value = %d\n", minimal_score[index_num]);
		}
		fprintf(bilan, "    Number of seeds = %d\n", seed_hits_gv);
		fprintf(bilan, "    Edges = %d", edges_gv);
		if (as_percent_gv)
			fprintf(bilan, " (as percent)\n");
		else
			fprintf(bilan, " (as integer)\n");
		fprintf(bilan, "    SW match = %ld\n", match);
		fprintf(bilan, "    SW mismatch = %ld\n", mismatch);
		fprintf(bilan, "    SW gap open penalty = %ld\n", gap_open);
		fprintf(bilan, "    SW gap extend penalty = %ld\n", gap_extension);
		fprintf(bilan, "    SW ambiguous nucleotide = %ld", score_N);
		if (score_N > 0) fprintf(bilan, " <-- %sWarning!%s Positive score set for ambiguous nucleotides.\n", "\033[0;33m", "\033[0m");
		else fprintf(bilan, "\n");
		if (yes_SQ)
			fprintf(bilan, "    SQ tags are output\n");
		else
			fprintf(bilan, "    SQ tags are not output\n");
#ifdef _OPENMP
		fprintf(bilan, "    Number of threads = %d\n", numcpu_gv);
#else
		fprintf(bilan, "    Number of threads = 1 (OpenMP is not supported with your current C++ compiler).\n");
#endif
		fprintf(bilan, "    Reads file = %s\n\n", inputreads);
		fclose(bilan);
	}

	// pointer to the split read
	// (the read which is split between any two file sections)
	char* split_read = NULL;
	// pointer to the position in the split read where to attach
	// the connecting part of the split read (and possibly its pair)
	char* split_read_ptr = NULL;
	// the number of lines to offset at the top of the current
	// file section
	int32_t offset_pair_from_top = 0;
	// map<reference sequence, vector<list of reads aligned to reference sequence> > otu_map
	map<string, vector<string>> otu_map;

	// Loop through all mmap'd read file sections
	while (file_s < file_sections)
	{
		if (samout_gv)
		{
			acceptedsam.open(acceptedstrings_sam, ios::app | ios::binary);
			if (!acceptedsam.good())
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not open SAM output file for writing: %s.\n",
					startColor, "\033[0m", __LINE__, __FILE__, acceptedstrings_sam);
				exit(EXIT_FAILURE);
			}
		}
		if (blastout_gv)
		{
			acceptedblast.open(acceptedstrings_blast, ios::app | ios::binary);
			if (!acceptedblast.good())
			{
				fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not open BLAST output file for writing: %s.\n",
					startColor, "\033[0m", __LINE__, __FILE__, acceptedstrings_blast);
				exit(EXIT_FAILURE);
			}
		}
		if (map_size_set_gv)
			eprintf("\n  %sBegin mmap reads section # %d%s:\n", "\033[4m", file_s + 1, "\033[0m");
		// begin file memory map
		TIME(s);
		char* raw = NULL;
		uint64_t strs = 0;
		char* finalnt = NULL;
		// the length of the split read in file part i+1 (from beginning of file part)
		uint32_t reads_offset_f = 0;
		// the length of the split read in file part i (from end of file part)
		uint32_t reads_offset_e = 0;
		char** reads = NULL;
		uint32_t min_lnwin = sizeof(uint32_t);
		for (int lwin = 0; lwin < num_databases; lwin++)
		{
			min_lnwin = (lnwin[lwin] < min_lnwin) ? lnwin[lwin] : min_lnwin;
		}
		if (map_size_set_gv)
		{
			reads = mmap_reads(partial_file_size,
				inputreads,
				offset_map,
				raw,
				filesig,
				file_s,
				file_sections,
				offset_pair_from_top,
				split_read_ptr,
				split_read,
				strs,
				finalnt,
				reads_offset_f,
				reads_offset_e,
				min_lnwin);
		}
		else
		{
			reads = load_reads(inputreads,
				raw,
				number_total_read,
				full_file_size,
				finalnt);
			strs = number_total_read * 2;
		}
		TIME(f);
		eprintf("  Time to load reads and set up pointers [%.2f sec]\n", (f - s));
		// array of bits to represent all reads
		// a bit representing an accepted read is set to 1
		vector<bool> read_hits(strs, false); // Read::hit
		// array of bits to represent all reads
		// a bit representing an accepted read with < %id 
		// and < %coverage is set to 0
		vector<bool> read_hits_denovo(strs, true); // Read::hit_denovo
		// array of uint16_t to represent all reads, if the read was aligned with a maximum SW score, its number of alignments is incremeted by 1
		uint16_t *read_max_SW_score = new uint16_t[strs]; // Read::max_SW_score
		memset(read_max_SW_score, 0, sizeof(uint16_t)*strs);
		// map accessed by read number, storing a pair <index for smallest SSW score, pointer to array of num_best_hits_gv>
		map<uint64_t, alignment_struct> read_hits_align_info; // Read::hits_align_info
		// number of alignments to output per read
		int32_t *num_alignments_x = NULL; // Read::
		// output num_alignments_gv alignments per read
		if (num_alignments_gv > 0)
		{
			num_alignments_x = new int32_t[strs];
			for (uint64_t s = 0; s < strs; s++) num_alignments_x[s] = num_alignments_gv;
		}

		// loop through every index passed to option --ref (ex. SSU 16S and SSU 18S)
		for (uint16_t index_num = 0; index_num < (uint16_t)myfiles.size(); index_num++)
		{
			// covert part number into a string
			stringstream prt_str;
			uint16_t part = 0;
			prt_str << part;
			string part_str = prt_str.str();
			eprintf("\n  Begin analysis of: %s%s%s\n", "\033[0;34m", (char*)(myfiles[index_num].first).c_str(), "\033[0m");
			if (file_s == 0)
			{
				eprintf("    Seed length = %d\n", lnwin[index_num]);
				eprintf("    Pass 1 = %d, Pass 2 = %d, Pass 3 = %d\n", skiplengths[index_num][0], skiplengths[index_num][1], skiplengths[index_num][2]);
				eprintf("    Gumbel lambda = %f\n", gumbel[index_num].first);
				eprintf("    Gumbel K = %f\n", gumbel[index_num].second);
				eprintf("    Minimal SW score based on E-value = %d\n", minimal_score[index_num]);
			}

			// for each partial file of burst trie index (part_0 .. part_x)
			for (part = 0; part < num_index_parts[index_num]; part++)
			{
				eprintf("    Loading index part %d/%u ... ", part + 1, num_index_parts[index_num]);
				TIME(s);
				// number of reference sequences to search alignments for before choosing the best one
				int32_t *best_x = NULL;
				// search min_lis_gv reference sequences for alignments
				if (min_lis_gv > 0)
				{
					best_x = new int32_t[strs];
					for (uint64_t s = 0; s < strs; s++)
						best_x[s] = min_lis_gv;
				}

				// memory buffer to store the reference sequence database
				char* buffer = NULL; // References
				// pointer to start of each sequence in buffer
				char** reference_seq = NULL; // References
				// length of each sequence in buffer
				uint64_t* reference_seq_len = NULL; // References

				// 9-mer look-up tables
				kmer *lookup_tbl = NULL; // Index
				// 19-mer position look-up tables
				kmer_origin* positions_tbl = NULL; // Index
				// number of elements in the table
				uint32_t number_elements = 0; // Index

				uint64_t seq_part_size = index_parts_stats_vec[index_num][part].seq_part_size; // References::load
				uint32_t numseq_part = index_parts_stats_vec[index_num][part].numseq_part; // References::load
				uint64_t start_part = index_parts_stats_vec[index_num][part].start_part; // References::load
#pragma omp master
				{
					// 2. load the index part (9-mer lookup table, mini-burst tries and positions table)
					load_index((char*)(myfiles[index_num].second).c_str(), part_str, lookup_tbl, positions_tbl, number_elements, lnwin[index_num]);

// References::load ------------------------------------------------------------------------------>
					// block of memory to hold all ids + reference sequences
					buffer = new char[(seq_part_size + 1)]();
					if (buffer == NULL)
					{
						fprintf(stderr, "    %sERROR%s: could not allocate memory for reference sequence buffer (paralleltraversal.cpp)\n", startColor, "\033[0m");
						exit(EXIT_FAILURE);
					}
					// pointer to the start of every sequence in the buffer
					reference_seq = new char*[(numseq_part << 1)]();
					if (reference_seq == NULL)
					{
						fprintf(stderr, "    %sERROR%s: Line %d: %s could not allocate memory for reference_seq\n",
							startColor, "\033[0m", __LINE__, __FILE__);
						exit(EXIT_FAILURE);
					}
					// length of every sequence in the buffer
					reference_seq_len = new uint64_t[numseq_part]();
					if (reference_seq_len == NULL)
					{
						fprintf(stderr, "    %sERROR%s: Line %d: %s could not allocate memory for reference_seq_len\n",
							startColor, "\033[0m", __LINE__, __FILE__);
						exit(EXIT_FAILURE);
					}
					// load the reference sequences for SW alignment
					load_ref((char*)(myfiles[index_num].first).c_str(),
						buffer,
						reference_seq,
						reference_seq_len,
						seq_part_size,
						numseq_part,
						start_part,
						1);
				}
// <------------------------------------------------------------ References::load

				TIME(f);
				eprintf(" done [%.2f sec]\n", (f - s));

				// begin the parallel traversal
				eprintf("    Begin index search ... ");
				TIME(s);

				uint32_t bit_vector_size = (partialwin[index_num] - 2) << 2; // on each index part. Init in Main, accessed in Worker
				uint32_t offset = (partialwin[index_num] - 3) << 2; // on each index part

				// only search the forward xor reverse strand
				int32_t max = 0;
				forward_gv = true;
				if (forward_gv ^ reverse_gv) max = 1;
				// search both strands
				else max = 2;

				// search the forward and/or reverse strands
				for (int32_t strand = 0; strand < max; strand++)
				{
					// loop through all of the reads in the file 
#pragma omp parallel for num_threads(numcpu_gv) shared(lookup_tbl,positions_tbl,buffer,reference_seq,reference_seq_len,read_hits_align_info,read_hits,read_max_SW_score) schedule(dynamic,256)
					for (uint64_t readn = 1; readn < strs; readn += 2)
					{
#ifdef debug_align
						cout << "readn = " << readn << endl; //TESTING
#endif

// parallelTraversalJob -------------------------------------------------------------------->
			// for reverse reads
						if (!forward_gv)
						{
							// output the first num_alignments_gv alignments
							if (num_alignments_gv > 0)
							{
								// all num_alignments_gv alignments have been output
								if (num_alignments_x[readn] < 0) continue;
							}
							// the maximum scoring alignment has been found, go to next read
							// (unless all alignments are being output)
							else if ((num_best_hits_gv > 0) &&
								(min_lis_gv > 0) &&
								(read_max_SW_score[readn] == num_best_hits_gv)) continue;
						}

						// read on integer alphabet {0,1,2,3}
						char myread[READLEN] = "";
						char* str = reads[readn];
						char* _myread = &myread[0];
						// keep record of what position on the read an ambiguous char exists,
						// the char at this position will be changed to '4' during
						// sequence alignment
						int32_t ambiguous_nt[READLEN] = { 0 };
						// size of ambiguous_nt
						uint32_t size_ambiguous_nt = 0;
						// flag to count only 1 alignment per read
						bool read_to_count = true;
						// length of read
						uint32_t readlen = 0;
						// change the read into an integer alphabet -- FASTA
						if (filesig == '>')
						{
							// str != finalnt (finalnt marks the end of mmap'd region)
							// str != '>' ('>' marks the start of the next read)
							// str != '\0' ('\0' marks the end of split-read, which
							// must be terminated by '\0')
							while ((str != finalnt) && (*str != '>') && (*str != '\0'))
							{
								if (*str != '\n' && *str != '\r')
								{
									*_myread = nt_table[(int)*str];
									if (*_myread == 4)
									{
										*_myread = 0;
										ambiguous_nt[size_ambiguous_nt++] = readlen;
									}
									_myread++;
									readlen++;
									if (readlen > READLEN)
									{
										fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] at least one of your reads is > %d nt \n",
											startColor, "\033[0m", __LINE__, __FILE__, READLEN);
										fprintf(stderr, "  Please check your reads or contact the authors.\n");
										exit(EXIT_FAILURE);
									}
								}
								str++;
							}
							*_myread = '\n';
						}
						// change the read into an integer alphabet -- FASTQ
						else
						{
							while (*str != '\n' && *str != '\r')
							{
								*_myread = nt_table[(int)*str++];
								if (*_myread == 4)
								{
									// searching the read using the universal levenshtein automaton
									// can only be done on 0-3 alphabet as burst trie is built 
									// on 0-3 alphabet
									*_myread = 0;
									ambiguous_nt[size_ambiguous_nt++] = readlen;
								}
								_myread++;
								readlen++;
								if (readlen > READLEN)
								{
									fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] at least one of your reads is > %d nt \n",
										startColor, "\033[0m", __LINE__, __FILE__, READLEN);
									fprintf(stderr, "  Please check your reads or contact the authors.\n");
									exit(EXIT_FAILURE);
								}
							}
							*_myread = '\n';
						}

						// find the minimum sequence length
						readlen < min_read_len ? min_read_len = readlen : min_read_len;
						// find the maximum sequence length
						readlen > max_read_len ? max_read_len = readlen : max_read_len;

						// the read length is too short
						if (readlen < lnwin[index_num])
						{
							fprintf(stderr, "\n  %sWARNING%s: At least one of the reads is shorter "
								"than %u nucleotides, ", "\033[0;33m", "\033[0m",
								lnwin[index_num]);
							fprintf(stderr, "by default it will not be searched\n ");
							continue;
						}

						// create the reverse strand
						if (!forward_gv)
						{
							// faster than xor algorithm
							char myread_rc[READLEN] = "";
							char* revcomp = &myread[readlen - 1];
							for (uint32_t j = 0; j < readlen; j++)
								myread_rc[j] = complement[(int)*revcomp--];
							myread_rc[readlen] = '\n';
							memcpy(&myread[0], &myread_rc[0], READLEN);
#ifdef debug_align
							cout << "read (before subst. 4's for N's): " << endl;
							char *cp = myread;
							while (*cp != '\n') cout << (int)*cp++;
							cout << endl;
#endif                      
						}//~if (REVERSE)

						// array of positions of window hits on the reference sequence
						vector< id_win > id_win_hits;
						// number of windows hit to the reference sequence(s)
						uint32_t readhit = 0;
						uint32_t windowshift = skiplengths[index_num][0];
						// keep track of windows which have been already traversed in the burst trie
						vector<bool> read_index_hits(readlen);
						// Pass number (possible value 0,1,2)
						uint32_t pass_n = 0;
						// the maximum SW score attainable for this read
						uint32_t max_SW_score = readlen*match;
						// loop for each new Pass to granulate seed search intervals
						bool search = true;
						std::vector<MYBITSET> vbitwindowsf;
						std::vector<MYBITSET> vbitwindowsr;
						do
						{
#ifdef debug_align
							cout << "\tpass = " << pass_n << endl; //TESTING
#endif          
							uint32_t numwin = (readlen - lnwin[index_num] + windowshift) / windowshift;
							uint32_t read_index = 0;
							// iterate over windows of the template string
							for (uint32_t win_num = 0; win_num < numwin; win_num++)
							{
								// skip position, seed at this position has already been searched for in a previous Pass
								if (read_index_hits[read_index]) goto check_score;
								// search position, set search bit to true
								else read_index_hits[read_index].flip();
								{
									// this flag it set to true if a match is found during
									// subsearch 1(a), to skip subsearch 1(b)
									bool accept_zero_kmer = false;
									// ids for k-mers that hit the database
									vector< id_win > id_hits;
									vbitwindowsf.resize(bit_vector_size);
									std::fill(vbitwindowsf.begin(), vbitwindowsf.end(), 0);
									//MYBITSET bitwindowsf[bit_vector_size]; // Bug 2. bit_vector_size is not known at compile time. Not supported in standard C++ => Fails in VS.
									//memset(&bitwindowsf[0],0,bit_vector_size);               
									// build the first bitvector window
									init_win_f(&myread[read_index + partialwin[index_num]],
										// [w_1] forward k = 1
										// bitwindows[0][0][0]
										&vbitwindowsf[0],
										// bitwindows[0][1][0]
										&vbitwindowsf[4],
										numbvs[index_num]);
									uint32_t keyf = 0;
									char *keyf_ptr = &myread[read_index];
									// build hash for first half windows (foward and reverse)
									for (uint32_t g = 0; g < partialwin[index_num]; g++)
										(keyf <<= 2) |= (uint32_t)*keyf_ptr++;
									// do traversal if the exact half window exists in the burst trie
									if ((lookup_tbl[keyf].count > minoccur) && (lookup_tbl[keyf].trie_F != NULL))
									{
										/* subsearch (1)(a) d([p_1],[w_1]) = 0 and d([p_2],[w_2]) <= 1;
										 *
										 *  w = |------ [w_1] ------|------ [w_2] ------|
										 *  p = |------ [p_1] ------|------ [p_2] ----| (0/1 deletion in [p_2])
										 *              or
										 *    = |------ [p_1] ------|------ [p_2] ------| (0/1 match/substitution in [p_2])
										 *        or
										 *    = |------ [p_1] ------|------ [p_2] --------| (0/1 insertion in [p_2])
										 *
										 */
#ifdef debug_align
										cout << "\tsearch forward mini-burst trie..\n"; //TESTING
#endif                   
										traversetrie_align(lookup_tbl[keyf].trie_F,
											0,
											0,
											// win2f_k1_ptr
											&vbitwindowsf[0],
											// win2f_k1_full
											&vbitwindowsf[offset],
											accept_zero_kmer,
											id_hits,
											readn,
											read_index,
											partialwin[index_num]);
#ifdef debug_align
										cout << "\tdone!\n"; //TESTING
#endif                      
									}//~if exact half window exists in the burst trie                 
									// only search if an exact match has not been found
									if (!accept_zero_kmer)
									{
										vbitwindowsr.resize(bit_vector_size);
										std::fill(vbitwindowsr.begin(), vbitwindowsr.end(), 0);
										// build the first bitvector window
										init_win_r(&myread[read_index + partialwin[index_num] - 1],
											&vbitwindowsr[0],
											&vbitwindowsr[4],
											numbvs[index_num]);
										uint32_t keyr = 0;
										char *keyr_ptr = &myread[read_index + partialwin[index_num]];
										// build hash for first half windows (foward and reverse)
										for (uint32_t g = 0; g < partialwin[index_num]; g++)
											(keyr <<= 2) |= (uint32_t)*keyr_ptr++;
										// continue subsearch (1)(b)
										if ((lookup_tbl[keyr].count > minoccur) && (lookup_tbl[keyr].trie_R != NULL))
										{
											/* subsearch (1)(b) d([p_1],[w_1]) = 1 and d([p_2],[w_2]) = 0;
											 *
											 *  w =    |------ [w_1] ------|------ [w_2] -------|
											 *  p =      |------- [p_1] ---|--------- [p_2] ----| (1 deletion in [p_1])
											 *              or
											 *    =    |------ [p_1] ------|------ [p_2] -------| (1 match/substitution in [p_1])
											 *        or
											 *    = |------- [p_1] --------|---- [p_2] ---------| (1 insertion in [p_1])
											 *
											 */
#ifdef debug_align
											cout << "\tsearch reverse mini-burst trie..\n"; //TESTING
#endif                     
											traversetrie_align(lookup_tbl[keyr].trie_R,
												0,
												0,
												/* win1r_k1_ptr */
												&vbitwindowsr[0],
												/* win1r_k1_full */
												&vbitwindowsr[offset],
												accept_zero_kmer,
												id_hits,
												readn,
												read_index,
												partialwin[index_num]);
#ifdef debug_align
											cout << "\tdone!\n"; //TESTING
#endif                        
										}//~if exact half window exists in the reverse burst trie                    
									}//~if (!accept_zero_kmer)                         
									// associate the ids with the read window number
									if (!id_hits.empty())
									{
										for (uint32_t i = 0; i < id_hits.size(); i++)
										{
											id_win_hits.push_back(id_hits[i]);
										}
										readhit++;
									}
								}
							check_score:
								// continue read analysis if threshold seeds were matched
								if (win_num == numwin - 1)
								{
									compute_lis_alignment(size_ambiguous_nt,
										readhit,
										id_win_hits,
										positions_tbl,
										read_max_SW_score,
										search,
										best_x,
										readn,
										num_alignments_x,
										readlen,
										lnwin[index_num],
										index_num,
										reference_seq_len,
										myread,
										ambiguous_nt,
										scoring_matrix,
										reference_seq,
										gap_open,
										gap_extension,
										minimal_score[index_num],
										read_hits,
										total_reads_mapped,
										reads_matched_per_db,
										part,
										read_hits_align_info,
										max_SW_score,
										read_to_count,
										total_reads_mapped_cov,
										read_hits_denovo,
										filesig,
										strs,
										file_s,
										file_sections,
										reads,
										finalnt,
										gumbel[index_num].first,
										gumbel[index_num].second,
										full_ref[index_num],
										full_read[index_num],
										acceptedblast,
										acceptedsam);
									// the read was not accepted at current window skip length,
									// decrease the window skip length
									if (search)
									{
										// last (3rd) Pass has been made
										if (pass_n == 2) search = false;
										else
										{
											// the next interval size equals to the current one, skip it
											while ((pass_n < 3) &&
												(skiplengths[index_num][pass_n] == skiplengths[index_num][pass_n + 1])) ++pass_n;
											if (++pass_n > 2) search = false;
											// set interval skip length for next Pass
											else windowshift = skiplengths[index_num][pass_n];
										}
									}
									// do not offset final window on read
									break;
								}//~( win_num == NUMWIN-1 )
								read_index += windowshift;
							}//~for (each window)                
						  //~while all three window skip lengths have not been tested, or a match has not been found
						} while (search);

						// the read didn't align (for --num_alignments [INT] option),
						// output null alignment string
						if (!read_hits[readn] && !forward_gv && (num_alignments_gv > -1))
						{
							// do not output read for de novo OTU clustering
							// (it did not pass the E-value threshold)
							if (de_novo_otu_gv && read_hits_denovo[readn]) read_hits_denovo[readn].flip();
							// output null alignment string
							if (print_all_reads_gv)
							{
#pragma omp critical
								{
									s_align* null_alignment = NULL;
									if (blastout_gv && blast_tabular)
									{
										report_blast(acceptedblast, // blast output file
											null_alignment, // SW alignment cigar
											reads[readn - 1] + 1, //read name
											0, // read sequence (in integer format)
											0, // read quality
											0, // reference name
											0, // reference sequence
											0, // e-value score
											0, // read length (to compute the masked regions)
											0, // bitscore
											0, // forward or reverse strand
											0, // %id
											0, // %query coverage
											0, // number of mismatches
											0); // number of gaps
									}
									if (samout_gv)
									{
										report_sam(acceptedsam, // sam output file
											null_alignment, // SW alignment cigar
											reads[readn - 1] + 1, // read name
											0, // read sequence (in integer format)
											0, // read quality
											0, // reference name
											0, // reference sequence
											0, // read length (to compute the masked regions)
											0, // forward or reverse strand
											0); // edit distance
									}
								}//~if print_all_reads_gv
							}// allow writing to file 1 thread at a time
						}//~if read didn't align
// <---------------------------------------------------------------------------- parallelTraversalJob
					}//~pragma omp for (each read)
#pragma omp barrier    
		  // search the reverse strand (default)
					forward_gv = false;
				}//for forward and/or reverse strands
				TIME(f);
				eprintf(" done [%.2f sec]\n", (f - s));
				eprintf("    Freeing index ... ");
				TIME(s);
#pragma omp master
				{
					// free the positions table
					if (positions_tbl != NULL)
					{
						for (uint32_t i = 0; i < number_elements; i++)
						{
							if (positions_tbl[i].arr != NULL)
							{
								delete[] positions_tbl[i].arr;
								positions_tbl[i].arr = NULL;
							}
						}
						delete[] positions_tbl;
						positions_tbl = NULL;
					}
					// free reference sequences loaded into memory
					if (buffer != NULL)
					{
						delete[] buffer;
						buffer = NULL;
					}
					if (reference_seq != NULL)
					{
						delete[] reference_seq;
						reference_seq = NULL;
					}
					if (reference_seq_len != NULL)
					{
						delete[] reference_seq_len;
						reference_seq_len = NULL;
					}
					// free 9-mer look-up tables and mini-burst tries
					if (lookup_tbl != NULL)
					{
						for (uint32_t i = 0; i < uint32_t(1 << lnwin[index_num]); i++)
						{
							if (lookup_tbl[i].trie_F != NULL)
							{
								delete[] lookup_tbl[i].trie_F;
								lookup_tbl[i].trie_F = NULL;
								lookup_tbl[i].trie_R = NULL;
							}
							if (lookup_tbl[i].trie_R != NULL)
							{
								delete[] lookup_tbl[i].trie_R;
								lookup_tbl[i].trie_R = NULL;
							}
						}
						delete[] lookup_tbl;
						lookup_tbl = NULL;
					}
				}
				TIME(f);
				eprintf(" done [%.2f sec]\n", (f - s));
				// increment the index part to next file
				prt_str.str("");
				prt_str << (part + 1);
				part_str = prt_str.str();
				// clear array holding number of reference sequences to
				// search per read prior to choosing best alignment
				if (best_x != NULL) delete[] best_x;
			}//~for all parts of the index            
		}//~for all indexes (provided as a list using option --ref)

		// clear array holding number of alignments output for each read
		if (num_alignments_x != NULL) delete[] num_alignments_x;
		delete[] read_max_SW_score;
		eprintf("    Total number of reads mapped (incl. all reads file sections searched): %llu\n", total_reads_mapped);

// writeAlignmentsJob -------------------------------------------------------------->
		// filter the sequences by %id and %query coverage, output them if --best INT
		if (min_lis_gv > -1)
		{
			if (samout_gv || blastout_gv) eprintf("    Writing alignments ... ");
#ifdef debug_mmap
			cout << "total index_num = " << myfiles.size() << endl;
			for (uint32_t index_num = 0; index_num < myfiles.size(); index_num++)
				cout << "\t parts = " << num_index_parts[index_num] << endl;
#endif            
			TIME(s);
			// loop through all databases
			for (uint32_t index_num = 0; index_num < myfiles.size(); index_num++)
			{
#ifdef debug_output
				cout << "index_num = " << index_num << endl;
#endif
				// loop through each section of a database
				for (uint16_t part = 0; part < num_index_parts[index_num]; part++)
				{
					uint64_t seq_part_size = index_parts_stats_vec[index_num][part].seq_part_size;
					uint32_t numseq_part = index_parts_stats_vec[index_num][part].numseq_part;
					uint64_t start_part = index_parts_stats_vec[index_num][part].start_part;
#ifdef debug_output
					cout << "part = " << part << endl;
					cout << "seq_part_size = " << seq_part_size << endl;
					cout << "numseq_part = " << numseq_part << endl;
					cout << "start_part = " << start_part << endl;
#endif                    
					// block of memory to hold all ids + reference sequences
					char* buffer = new char[(seq_part_size + 1)]();
					if (buffer == NULL)
					{
						fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for reference sequence buffer\n",
							startColor, "\033[0m", __LINE__, __FILE__);
						exit(EXIT_FAILURE);
					}

					// pointer to the start of every sequence in the buffer
					char** reference_seq = new char*[(numseq_part << 1)]();
					if (reference_seq == NULL)
					{
						fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for reference_seq\n",
							startColor, "\033[0m", __LINE__, __FILE__);
						exit(EXIT_FAILURE);
					}

					// length of every sequence in the buffer (is not required here)
					uint64_t* reference_seq_len = NULL;
					load_ref((char*)(myfiles[index_num].first).c_str(),
						buffer,
						reference_seq,
						reference_seq_len,
						seq_part_size,
						numseq_part,
						start_part,
						0);

					// run through all the reads, output those which aligned
					for (uint64_t readn = 1; readn < strs; readn += 2)
					{
						map<uint64_t, alignment_struct >::iterator alignment = read_hits_align_info.find(readn);
						// this read does not have any alignment
						if (alignment == read_hits_align_info.end())
						{
							// (output NULL alignment only once)
							if ((index_num == 0) && (part == 0))
							{
								// do not output this read for de novo clustering
								// (it did not pass the E-value threshold)
								if (de_novo_otu_gv && read_hits_denovo[readn]) read_hits_denovo[readn].flip();
								// output null string for read alignment
								if (print_all_reads_gv)
								{
									s_align* null_alignment = NULL;
									if (blastout_gv && blast_tabular)
									{
										report_blast(acceptedblast, // blast output file
											null_alignment, // SW alignment cigar
											reads[readn - 1] + 1, //read name
											0, // read sequence (in integer format)
											0, // read quality
											0, // reference name
											0, // reference sequence
											0, // e-value score
											0, // read length (to compute the masked regions)
											0, // bitscore
											0, // forward or reverse strand
											0, // %id
											0, // %query coverage
											0, // number of mismatches
											0); // number of gaps
									}
									if (samout_gv)
									{
										report_sam(acceptedsam, // sam output file
											null_alignment, // SW alignment cigar
											reads[readn - 1] + 1, // read name
											0, // read sequence (in integer format)
											0, // read quality
											0, // reference name
											0, // reference sequence
											0, // read length (to compute the masked regions)
											0, // forward or reverse strand
											0); // edit distance
									}
								}
							}
							// go to next read
							continue;
						}
#ifdef debug_output
						cout << "readn = " << readn << endl;
#endif
						// get all alignments for this read
						s_align* ptr_alignment = alignment->second.ptr;
						if (ptr_alignment == NULL)
						{
							fprintf(stderr, "  ERROR: s_align* ptr_alignment == NULL, this should not be possible.\n");
							exit(EXIT_FAILURE);
						}
#ifdef debug_output
						cout << "ptr_alignment->index_num = " << (uint16_t)ptr_alignment->index_num << endl;
						cout << "ptr_alignment->score1 = " << (int32_t)ptr_alignment->score1 << endl;
						cout << "ptr_alignment->part = " << (uint16_t)ptr_alignment->part << endl;
#endif
						// OTU-map: index of alignment holding maximum SW score
						uint32_t index_max_score = alignment->second.max_index;
						// loop through all of the best alignments for this read
						for (uint32_t p = 0; p < alignment->second.size; p++)
						{
#ifdef debug_output
							cout << "best_hit = " << p << endl;
#endif
							// continue loop if the reference sequence in this alignment
							// belongs to the database section currently loaded into RAM
							if ((ptr_alignment->index_num == index_num) && (ptr_alignment->part == part))
							{
								// format read & get read length
								char myread[READLEN];
								uint32_t readlen = ptr_alignment->readlen;
								// format forward read from char to int
								if (ptr_alignment->strand)
								{
									format_forward(reads[readn], &myread[0], filesig);
								}
								// format reverse-complement read from char to int
								else
								{
									char* end_read = NULL;
									// FASTA
									if (filesig == '>')
									{
										// if split-read (for file_s > 0)
										// note: the split-read contains two reads (a pair, in case reads are paired) that
										//   are located at another address space than the contiguous memory mapped
										//   section of reads
										if ((readn < 4) && (file_s > 0))
										{
#ifdef debug_output
											cout << "process split-read" << endl; //TESTING
#endif
											end_read = reads[readn];
											while (*end_read != '\0') end_read++;
											if (*end_read == '\0') end_read--; //account for '\0'
											if (*end_read == '\n' || *end_read == '\r') end_read--; //account for '\n'
										}
										// last read in the partial file section
										else if (readn >= (strs - 2))
										{
											end_read = reads[readn];
											// if processing last file section, the final read will end with '\0'
											if (file_s == file_sections - 1)
											{
#ifdef debug_output
												cout << "process last read in final file section .." << endl; //TESTING
#endif
												while (*end_read++ != '\0');
												// back-track 3 nucleotides (\n\0_)
												if (*(end_read - 2) == '\n' || *(end_read - 2) == '\r') end_read--;
												// back-track 2 nucleotides (\0_)
												end_read -= 2;
											}
											// if processing a file section > 0 and < last file section, the final read will end with '>' (beginning of split-read)
											else
											{
#ifdef debug_output
												cout << "process last read in file section > 0 and < final file section .." << endl; //TESTING
#endif
												while (*end_read++ != '>');
												// back-track 3 nucleotides (\n>_)
												end_read -= 3;
											}
										}
										// all other reads
										else end_read = reads[readn + 1] - 2;
									}
									// FASTQ
									else end_read = reads[readn] + readlen - 1;
									format_rev(reads[readn], end_read, &myread[0], filesig);
								}//~reverse-complement read
								// get the edit distance between reference and read
								uint32_t ref_seq = ptr_alignment->ref_seq;
								double id = 0;
								char to_char[5] = { 'A','C','G','T','N' };
								char* ref_seq_ptr = reference_seq[(2 * ref_seq) + 1];
								char* read_seq_ptr = myread;
								int32_t qb = ptr_alignment->ref_begin1;
								int32_t pb = ptr_alignment->read_begin1;
								uint32_t mismatches = 0;
								uint32_t gaps = 0;
#ifdef debug_output
								cout << "ref_seq = " << ref_seq << endl;
								cout << "ptr_alignment->ref_begin1 = " << ptr_alignment->ref_begin1 << endl;
								cout << "ptr_alignment->read_begin1 = " << ptr_alignment->read_begin1 << endl;
#endif                           
								for (uint32_t c2 = 0; c2 < ptr_alignment->cigarLen; ++c2)
								{
									uint32_t letter = 0xf & *(ptr_alignment->cigar + c2);
									uint32_t length = (0xfffffff0 & *(ptr_alignment->cigar + c2)) >> 4;
									if (letter == 0)
									{
										for (uint32_t u = 0; u < length; ++u)
										{
											if ((char)to_char[(int)*(ref_seq_ptr + qb)] != (char)to_char[(int)*(read_seq_ptr + pb)]) ++mismatches;
											else ++id;
											++qb;
											++pb;
										}
									}
									else if (letter == 1)
									{
										pb += length;
										gaps += length;
									}
									else
									{
										qb += length;
										gaps += length;
									}
								}
								int32_t align_len = abs(ptr_alignment->read_end1 + 1 - ptr_alignment->read_begin1);
								int32_t total_pos = mismatches + gaps + id;
								stringstream ss;
								ss.precision(3);
								ss << (double)id / total_pos << ' ' << (double)align_len / readlen;
								double align_id_round = 0.0;
								double align_cov_round = 0.0;
								ss >> align_id_round >> align_cov_round;
								//#define debug_id_cov
#ifdef debug_id_cov
								cout << "read tag: ";
								char* tt = reads[readn - 1];
								while (*tt != '\n') cout << (char)*tt++;
								cout << endl;
								cout << "align_len = " << align_len << endl;
								cout << "total_pos = " << total_pos << endl;
								cout << "id = " << id << endl;
								cout << "Score = " << ptr_alignment->score1 << endl;
								cout << "%id = " << (double)id / total_pos << endl;
								cout << "%cov = " << (double)align_len / readlen << endl;
								cout << "align_cov = " << (double)align_cov << endl;
								cout << "align_id = " << (double)align_id << endl;
								cout << "align_id_round = " << (double)align_id_round << endl;
								cout << "align_cov_round = " << (double)align_cov_round << endl;
#endif                              
								// alignment with the highest SW score passed
								// %id and %coverage thresholds
								if ((p == index_max_score) &&
									(align_id_round >= align_id) &&
									(align_cov_round >= align_cov))
								{
									// increment number of reads passing identity
									// and coverage threshold
									total_reads_mapped_cov++;
									// do not output read for de novo OTU construction
									// (it passed the %id/coverage thresholds)
									if (de_novo_otu_gv && read_hits_denovo[readn]) read_hits_denovo[readn].flip();
									// fill OTU map with highest-scoring alignment for the read
									if (otumapout_gv)
									{
										// reference sequence identifier for mapped read
										char ref_seq_arr[4000] = "";
										char* ref_seq_arr_ptr = ref_seq_arr;
										char* ref_seq_id_ptr = reference_seq[(2 * ref_seq)] + 1;
										while ((*ref_seq_id_ptr != ' ') && (*ref_seq_id_ptr != '\n') && (*ref_seq_id_ptr != '\r')) *ref_seq_arr_ptr++ = *ref_seq_id_ptr++;
										string ref_seq_str = ref_seq_arr;
										// read identifier
										char read_seq_arr[4000] = "";
										char* read_seq_arr_ptr = read_seq_arr;
										char* read_seq_id_ptr = reads[readn - 1] + 1;
										while ((*read_seq_id_ptr != ' ') && (*read_seq_id_ptr != '\n') && (*read_seq_id_ptr != '\r')) *read_seq_arr_ptr++ = *read_seq_id_ptr++;
										string read_seq_str = read_seq_arr;
										otu_map[ref_seq_str].push_back(read_seq_str);
									}
								}

								// output alignment to SAM or Blast-like formats
								if (samout_gv || blastout_gv)
								{
									// quality for FASTQ
									char* read_qual = NULL;
									if (filesig == '@')
									{
										// forward
										if (ptr_alignment->strand)
										{
											// if second part of split read or last read in file
											if (((readn == 3) && (file_s > 0)) || (readn >= (strs - 2)))
											{
#ifdef debug_output
												if (readn >= (strs - 2)) cout << "get quality for last (forward) read in file section\n"; //TESTING
#endif
												read_qual = reads[readn];
												int8_t numnewlines = 0;
												while (numnewlines < 2) { if (*read_qual++ == '\n' || *read_qual++ == '\r') numnewlines++; }
											}
											else read_qual = reads[readn + 1] - readlen - 1;
										}
										// reverse-complement
										else
										{
											// if second part of split read or last read in file
											if (((readn == 3) && (file_s > 0)) || (readn >= (strs - 2)))
											{
												read_qual = reads[readn];

												// last file section
												if (file_s == file_sections - 1)
												{
#ifdef debug_output
													if (readn >= (strs - 2)) cout << "get quality for last (reverse) read in last file section\n"; //TESTING
#endif
													while (*read_qual != '\0') read_qual++;
													if (*(read_qual - 3) == '\n' || *(read_qual - 3) == '\r') read_qual--;
													// account for '\n\0'
													read_qual -= 2;
												}
												// file section > 0 and < last file section
												else
												{
#ifdef debug_output
													if (readn >= (strs - 2)) cout << "get quality for last (reverse) read in (not last) file section\n"; //TESTING
#endif
													while (read_qual != finalnt) read_qual++;
													read_qual--;
												}
											}
											else read_qual = reads[readn + 1] - 2;
										}
									}//~if filesig == '@'

									if (blastout_gv)
									{
										uint32_t bitscore = (uint32_t)((float)((gumbel[index_num].first)*(ptr_alignment->score1) - log(gumbel[index_num].second)) / (float)log(2));
										double evalue_score = (double)(gumbel[index_num].second)*full_ref[index_num] * full_read[index_num] * pow(EXP, (-(gumbel[index_num].first)*ptr_alignment->score1));
										report_blast(acceptedblast, //blast output file
											ptr_alignment, //SW alignment cigar
											reads[readn - 1] + 1, //read name
											myread, //read sequence (in integer format)
											read_qual, //read quality
											reference_seq[(2 * ref_seq)] + 1, //reference name
											reference_seq[(2 * ref_seq) + 1], //reference sequence
											evalue_score, //e-value score
											readlen, //read length (to compute the masked regions)
											bitscore,
											ptr_alignment->strand,
											(double)s / total_pos, // %id
											(double)align_len / readlen, // %query coverage
											mismatches,
											gaps
										);
									}
									if (samout_gv)
									{
										report_sam(acceptedsam, //sam output file
											ptr_alignment, //SW alignment cigar
											reads[readn - 1] + 1, //read name
											myread,//read sequence (in integer format)
											read_qual, //read quality
											reference_seq[(2 * ref_seq)] + 1, //reference name
											reference_seq[(2 * ref_seq) + 1], //reference sequence
											readlen, //read length (to compute the masked regions)
											ptr_alignment->strand,
											mismatches + gaps); //edit distance
									}
								}//~if (samout_gv || blastout_gv)
							}//~if alignment at current database and index part loaded in RAM
							ptr_alignment++;
						}//~for all best alignments
					}//~for all reads

					// free buffer
					if (buffer != NULL)
					{
						delete[] buffer;
						buffer = NULL;
					}
					// free pointers to reference sequences
					if (reference_seq != NULL)
					{
						delete[] reference_seq;
						reference_seq = NULL;
					}
				}//~for every database section
			}//~for every database

			// free alignment information for all aligned reads
			for (uint64_t readn = 1; readn < strs; readn += 2)
			{
				map<uint64_t, alignment_struct >::iterator alignment = read_hits_align_info.find(readn);
				// this read does not have any alignment
				if (alignment != read_hits_align_info.end())
				{
					s_align* ptr_alignment = alignment->second.ptr;
					for (uint32_t p = 0; p < alignment->second.size; p++)
					{
						free(ptr_alignment->cigar);
						ptr_alignment->cigar = NULL;
						if (p + 1 < (uint32_t)num_best_hits_gv)
						{
							ptr_alignment++;
							// check whether an alignment exists
							if (ptr_alignment->cigar == NULL) break;
						}
					}
					// free memory for all alignments of this read
					delete[] alignment->second.ptr;
					alignment->second.ptr = NULL;
				}
			}
			TIME(f);
			if (samout_gv || blastout_gv) eprintf(" done [%.2f sec]\n", (f - s));
		}// if (min_lis_gv > -1)

		if (align_cov || align_id)
		{
			eprintf("    Total number of reads mapped with");
			if (align_id > 0) eprintf(" >= %.2lf identity,", align_id);
			if (align_cov > 0) eprintf(" >= %.2lf query coverage", align_cov);
			eprintf(" (incl. all reads file sections searched): %llu\n", total_reads_mapped_cov);
		}
// <------------------------------------------------------------------- writeAlignmentsJob

		if (blastout_gv)
		{
			if (acceptedblast.is_open()) acceptedblast.close();
			else
			{
				fprintf(stderr, "  %sERROR%s: file %s was not opened for writing.\n", startColor, "\033[0m", acceptedstrings_blast);
				exit(EXIT_FAILURE);
			}
		}

		if (samout_gv)
		{
			if (acceptedsam.is_open()) acceptedsam.close();
			else
			{
				fprintf(stderr, "  %sERROR%s: file %s was not opened for writing.\n", startColor, "\033[0m", acceptedstrings_sam);
				exit(EXIT_FAILURE);
			}
		}

		// output aligned and non-aligned reads to FASTA/FASTQ file
		report_fasta(acceptedstrings,
			ptr_filetype_or,
			ptr_filetype_ar,
			reads,
			strs,
			read_hits,
			file_s,
			finalnt);

		// output aligned and non-aligned reads with < %id and
		// < %coverage to FASTA/FASTQ file for de novo analysis
		if (de_novo_otu_gv)
		{
			// count number of reads output for de novo clustering
			for (uint64_t d = 1; d < strs; d += 2)
			{
				if (read_hits_denovo[d])
					total_reads_denovo_clustering++;
			}

			report_denovo(denovo_otus_file,
				reads,
				strs,
				read_hits_denovo,
				file_s,
				finalnt);
		}
		read_hits.clear();
		read_hits_denovo.clear();
		// free the split_read
		if (split_read != NULL)
		{
			delete[] split_read;
			split_read = NULL;
			split_read_ptr = NULL;
		}
		// record the start of the split_read if it exists
		if (file_s < file_sections - 1)
		{
			split_read = new char[(READLEN * 2)];
			if (split_read == NULL)
			{
				fprintf(stderr, "  %sERROR%s: could not allocate memory for the bridged read\n", startColor, "\033[0m");
				exit(EXIT_FAILURE);
			}
			// compute the first half of the split_read
			char *start = &raw[partial_file_size] - reads_offset_e - 1;
			char *end = &raw[partial_file_size];
			split_read_ptr = split_read;
#ifdef debug_mmap
			cout << "split read start: "; //TESTING
#endif     
			while (start != end)
			{
#ifdef debug_mmap
				cout << (char)*start; //TESTING
#endif
				*split_read_ptr++ = *start++;
			}
#ifdef debug_mmap
			cout << ".STOP." << endl; //TESTING
#endif  
		}// (s < file_sections - 1)
		if (map_size_set_gv)
		{
			// free the mmap'd file section
			unmmap_reads(raw, partial_file_size);
			offset_map += map_size_gv;
		}
		else
		{
			delete[] raw;
		}
		// last section of the full file, count strings until EOF is reached
		if (++file_s == file_sections - 1)
			partial_file_size = last_part_size;
		delete[] reads;
		reads = NULL;
	}//~while ( file_s < file_sections )

	// output OTU map to file
	if (otumapout_gv)
	{
		ofstream outfile(acceptedotumap_file, ios::app | ios::binary);
		map<string, vector<string> >::iterator otu_map_it;
		for (otu_map_it = otu_map.begin(); otu_map_it != otu_map.end(); otu_map_it++)
		{
			// output the ref ID
			outfile << otu_map_it->first;
			// output all reads mapping to ref ID
			for (uint32_t i = 0; i < otu_map_it->second.size(); i++) outfile << "\t" << otu_map_it->second[i];
			outfile << "\n";
		}
		if (outfile.is_open()) outfile.close();
		else
		{
			fprintf(stderr, "  %sERROR%s: file %s was not opened for writing.\n", startColor, "\033[0m", acceptedotumap_file);
			exit(EXIT_FAILURE);
		}
		// free memory for OTU mapping file
		if (acceptedotumap_file != NULL)
		{
			delete[] acceptedotumap_file;
			acceptedotumap_file = NULL;
		}
	}

	// free scoring matrix
	free(scoring_matrix);
	scoring_matrix = NULL;

	// create a bilan (log file)
	if ((ptr_filetype_ar != NULL) && logout_gv)
	{
		FILE* bilan = fopen(logoutfile, "ab");
		if (bilan == NULL)
		{
			fprintf(stderr, "  %sERROR%s: could not open file %s \n", startColor, "\033[0m", logoutfile);
			exit(EXIT_FAILURE);
		}
		// output total number of reads
		fprintf(bilan, " Results:\n");
		fprintf(bilan, "    Total reads = %llu\n", number_total_read);
		if (de_novo_otu_gv)
		{
			fprintf(bilan, "    Total reads for de novo clustering = %llu\n", total_reads_denovo_clustering);
		}
		// output total non-rrna + rrna reads
		fprintf(bilan, "    Total reads passing E-value threshold = %llu (%.2f%%)\n", total_reads_mapped, (float)((float)total_reads_mapped / (float)number_total_read) * 100);
		fprintf(bilan, "    Total reads failing E-value threshold = %llu (%.2f%%)\n", number_total_read - total_reads_mapped, (1 - ((float)((float)total_reads_mapped / (float)number_total_read))) * 100);
		fprintf(bilan, "    Minimum read length = %u\n", min_read_len);
		fprintf(bilan, "    Maximum read length = %u\n", max_read_len);
		fprintf(bilan, "    Mean read length = %u\n", mean_read_len);
		fprintf(bilan, " By database:\n");
		// output stats by database
		for (uint32_t index_num = 0; index_num < myfiles.size(); index_num++)
		{
			fprintf(bilan, "    %s\t\t%.2f%%\n", (char*)(myfiles[index_num].first).c_str(), (float)((float)reads_matched_per_db[index_num] / (float)number_total_read) * 100);
		}
		if (otumapout_gv)
		{
			fprintf(bilan, " Total reads passing %%id and %%coverage thresholds = %llu\n", total_reads_mapped_cov);
			fprintf(bilan, " Total OTUs = %lu\n", otu_map.size());
		}
		time_t q = time(0);
		struct tm * now = localtime(&q);
		fprintf(bilan, "\n %s\n", asctime(now));
		fclose(bilan);
		// free memory of accepted strings
		if (logoutfile != NULL)
		{
			delete[] logoutfile;
			logoutfile = NULL;
		}
	}
	else if (otumapout_gv) otu_map.clear();

	// free memory of accepted strings
	if (acceptedstrings != NULL)
	{
		delete[] acceptedstrings;
		acceptedstrings = NULL;
	}
	if (acceptedstrings_sam != NULL)
	{
		delete[] acceptedstrings_sam;
		acceptedstrings_sam = NULL;
	}
	if (acceptedstrings_blast != NULL)
	{
		delete[] acceptedstrings_blast;
		acceptedstrings_blast = NULL;
	}
	if (denovo_otus_file != NULL)
	{
		delete[] denovo_otus_file;
		denovo_otus_file = NULL;
	}
	return;
}//~paralleltraversal()

void writeLog(Runopts & opts, Index & index, Readstats & readstats, Output & output)
{
	//FILE* bilan = fopen(output.logoutfile, "ab"); // output::logoutfile
	output.logout.open(output.logoutfile, std::ofstream::binary | std::ofstream::app);

	// output total number of reads
	output.logout << " Results:\n";
	output.logout << "    Total reads = " << readstats.number_total_read << "\n";
	//fprintf(bilan, " Results:\n");
	//fprintf(bilan, "    Total reads = %llu\n", readstats.number_total_read); // Readstats::number_total_read
	if (de_novo_otu_gv)
	{
		// total_reads_denovo_clustering = sum of all reads that have read::hit_denovo == true
		// either query DB or store in Readstats::total_reads_denovo_clustering
		output.logout << "    Total reads for de novo clustering = " << readstats.total_reads_denovo_clustering << "\n";
		//fprintf(bilan, "    Total reads for de novo clustering = %llu\n", total_reads_denovo_clustering);
	}
	// output total non-rrna + rrna reads
	output.logout << std::setprecision(2) << std::fixed;
	output.logout << "    Total reads passing E-value threshold = " << readstats.total_reads_mapped 
		<< " (" << (float)((float)readstats.total_reads_mapped / (float)readstats.number_total_read) * 100 << ")\n";
	output.logout << "    Total reads failing E-value threshold = "
		<< readstats.number_total_read - readstats.total_reads_mapped
		<< " ("	<< (1 - ((float)((float)readstats.total_reads_mapped / (float)readstats.number_total_read))) * 100 << ")\n";
	output.logout << "    Minimum read length = " << readstats.min_read_len << "\n";
	output.logout << "    Maximum read length = " << readstats.max_read_len << "\n";
	output.logout << "    Mean read length    = " << readstats.full_read_main / readstats.number_total_read << "\n";

	output.logout << " By database:\n";
	//fprintf(bilan, "    Total reads passing E-value threshold = %llu (%.2f%%)\n", 
	//	readstats.total_reads_mapped, (float)((float)readstats.total_reads_mapped / (float)readstats.number_total_read) * 100);
	//fprintf(bilan, "    Total reads failing E-value threshold = %llu (%.2f%%)\n", 
	//	readstats.number_total_read - readstats.total_reads_mapped, 
	//	(1 - ((float)((float)readstats.total_reads_mapped / (float)readstats.number_total_read))) * 100);
	//fprintf(bilan, "    Minimum read length = %u\n", readstats.min_read_len);
	//fprintf(bilan, "    Maximum read length = %u\n", readstats.max_read_len);
	//fprintf(bilan, "    Mean read length = %u\n", readstats.full_read_main / readstats.number_total_read);
	//fprintf(bilan, " By database:\n");
	// output stats by database
	for (uint32_t index_num = 0; index_num < opts.indexfiles.size(); index_num++)
	{
		output.logout << "    " << opts.indexfiles[index_num].first << "\t\t"
			<< (float)((float)readstats.reads_matched_per_db[index_num] / (float)readstats.number_total_read) * 100 << "\n";
		//fprintf(bilan, "    %s\t\t%.2f%%\n", (char*)(opts.indexfiles[index_num].first).c_str(), 
		//	(float)((float)reads_matched_per_db[index_num] / (float)readstats.number_total_read) * 100);
	}

	if (otumapout_gv)
	{
		output.logout << " Total reads passing %%id and %%coverage thresholds = " << readstats.total_reads_mapped_cov << "\n";
		output.logout << " Total OTUs = " << readstats.otu_total << "\n"; // otu_map.size()
		//fprintf(bilan, " Total reads passing %%id and %%coverage thresholds = %llu\n", total_reads_mapped_cov);
		//fprintf(bilan, " Total OTUs = %lu\n", otu_map.size());
	}
	time_t q = time(0);
	struct tm * now = localtime(&q);
	output.logout << "\n " << asctime(now) << "\n";
	output.logout.close();
} // ~writeLog

void writeAlignments(Runopts & opts, Index & index, References & refs, Readstats & readstats)
{

	if (samout_gv || blastout_gv)
		std::cout << "    Writing alignments ... \n";
          
	//TIME(s);
	// loop through all databases
	for (uint32_t index_num = 0; index_num < opts.indexfiles.size(); index_num++)
	{
		// loop through each section of a database
		for (uint16_t part = 0; part < index.num_index_parts[index_num]; part++)
		{
			//uint64_t seq_part_size = index.index_parts_stats_vec[index.index_num][index.part].seq_part_size;
			//uint32_t numseq_part = index.index_parts_stats_vec[index.index_num][index.part].numseq_part;
			//uint64_t start_part = index.index_parts_stats_vec[index.index_num][index.part].start_part;
                
			// block of memory to hold all ids + reference sequences
			//char* buffer = new char[(seq_part_size + 1)](); // References::buffer  index.index_parts_stats_vec[index.index_num][index.part].seq_part_size
			//if (buffer == NULL)
			//{
			//	fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for reference sequence buffer\n",
			//		startColor, "\033[0m", __LINE__, __FILE__);
			//	exit(EXIT_FAILURE);
			//}

			// pointer to the start of every sequence in the buffer
			//char** reference_seq = new char*[(numseq_part << 1)](); // index.index_parts_stats_vec[index.index_num][index.part].numseq_part;
			//if (reference_seq == NULL)
			//{
			//	fprintf(stderr, "  %sERROR%s: [Line %d: %s] could not allocate memory for reference_seq\n",
			//		startColor, "\033[0m", __LINE__, __FILE__);
			//	exit(EXIT_FAILURE);
			//}

			// length of every sequence in the buffer (is not required here)
			//uint64_t* reference_seq_len = NULL;
			//load_ref((char*)(myfiles[index_num].first).c_str(),
			//	buffer,
			//	reference_seq,
			//	reference_seq_len,
			//	seq_part_size,
			//	numseq_part,
			//	start_part,
			//	0);

			// run through all the reads, output those which aligned
			//for (uint64_t readn = 1; readn < strs; readn += 2)
			//writeAlignmentJob
			for (int readn = 0; readn < readstats.number_total_read; ++readn) {} //~for all reads

			 // free buffer
			//if (buffer != NULL)
			//{
			//	delete[] buffer;
			//	buffer = NULL;
			//}
			// free pointers to reference sequences
			//if (reference_seq != NULL)
			//{
			//	delete[] reference_seq;
			//	reference_seq = NULL;
			//}
		}//~for every database section
	}//~for every database

	 // free alignment information for all aligned reads
	//for (uint64_t readn = 1; readn < strs; readn += 2)
	//{
	//	map<uint64_t, alignment_struct >::iterator alignment = read_hits_align_info.find(readn);
		// this read does not have any alignment
	//	if (alignment != read_hits_align_info.end())
	//	{
	//		s_align* ptr_alignment = alignment->second.ptr;
	//		for (uint32_t p = 0; p < alignment->second.size; p++)
	//		{
	//			free(ptr_alignment->cigar);
	//			ptr_alignment->cigar = NULL;
	//			if (p + 1 < (uint32_t)num_best_hits_gv)
	//			{
	//				ptr_alignment++;
					// check whether an alignment exists
	//				if (ptr_alignment->cigar == NULL) break;
	//			}
	//		}
			// free memory for all alignments of this read
	//		delete[] alignment->second.ptr;
	//		alignment->second.ptr = NULL;
	//	}
	//}
	//TIME(f);
	if (samout_gv || blastout_gv) 
		std::cout << " Done ... \n";
		//eprintf(" done [%.2f sec]\n", (f - s));
} // ~writeAlignments

void writeAlignmentJob(
	Index & index, 
	References & refs, 
	Output & output, 
	Readstats & readstats, 
	Read & read
	)
{
	// TODO: remove this part. [jenya's message of 20171121 8:47]
	//map<uint64_t, alignment_struct >::iterator alignment = read_hits_align_info.find(readn); // 
	// this read does not have any alignment
	//if (alignment == read_hits_align_info.end())
	//if (read.hits_align_info.alignv.size() == 0)
	//{
		// (output NULL alignment only once)
	//	if ((index.index_num == 0) && (index.part == 0))
	//	{
			// do not output this read for de novo clustering
			// (it did not pass the E-value threshold)
	//		if (de_novo_otu_gv && read.hit_denovo) read.hit_denovo != read.hit_denovo; // flip read_hits_denovo[readn]
			// output null string for read alignment
	//		if (print_all_reads_gv)
	//		{
	//			s_align* null_alignment = NULL;
	//			if (blastout_gv && blast_tabular)
	//			{
	//				report_blast(
	//					output.acceptedblast, // blast output file
	//					null_alignment, // SW alignment cigar
	//					read.header.c_str(), //read name
	//					0, // read sequence (in integer format)
	//					0, // read quality
	//					0, // reference name
	//					0, // reference sequence
	//					0, // e-value score
	//					0, // read length (to compute the masked regions)
	//					0, // bitscore
	//					0, // forward or reverse strand
	//					0, // %id
	//					0, // %query coverage
	//					0, // number of mismatches
	//					0 // number of gaps
	//				);
	//			}
	//			if (samout_gv)
	//			{
	//				report_sam(
	//					output.acceptedsam, // sam output file
	//					null_alignment, // SW alignment cigar
	//					read.header.c_str(), // read name
	//					0, // read sequence (in integer format)
	//					0, // read quality
	//					0, // reference name
	//					0, // reference sequence
	//					0, // read length (to compute the masked regions)
	//					0, // forward or reverse strand
	//					0 // edit distance
	//				);
	//			}
	//		}
	//	}
		// go to next read
		//continue;
	//	return;
	//} // ~ if read has no alignment information

	// get all alignments for this read
	//s_align* ptr_alignment = alignment->second.ptr; // &read.hits_align_info.alignv[0];
	//if (ptr_alignment == NULL)
	//{
	//	fprintf(stderr, "  ERROR: s_align* ptr_alignment == NULL, this should not be possible.\n");
	//	exit(EXIT_FAILURE);
	//}

	// OTU-map: index of alignment holding maximum SW score
	uint32_t index_max_score = read.hits_align_info.max_index;
	// loop through all of the best alignments for this read
	//for (uint32_t p = 0; p < alignment->second.size; p++)
	for (uint32_t p = 0; p < read.hits_align_info.alignv.size(); ++p)
	{
		// continue loop if the reference sequence in this alignment
		// belongs to the database section currently loaded into RAM
		if ((read.hits_align_info.alignv[p].index_num == index.index_num) && (read.hits_align_info.alignv[p].part == index.part))
		{
			// format read & get read length
			//char myread[READLEN];
			uint32_t readlen = read.hits_align_info.alignv[p].readlen;
			// format forward read from char to int
			//if (read.hits_align_info.alignv[p].strand)
			//{
			//	format_forward(reads[readn], &myread[0], filesig); // read.isequence
			//}
			// format reverse-complement read from char to int
			//else
			//{
			//	char* end_read = NULL;
				// FASTA
			//	if (filesig == '>')
			//	{
					// if split-read (for file_s > 0)
					// note: the split-read contains two reads (a pair, in case reads are paired) that
					//   are located at another address space than the contiguous memory mapped
					//   section of reads
			//		if ((readn < 4) && (file_s > 0))
			//		{
			//			end_read = reads[readn];
			//			while (*end_read != '\0') end_read++;
			//			if (*end_read == '\0') end_read--; //account for '\0'
			//			if (*end_read == '\n' || *end_read == '\r') end_read--; //account for '\n'
			//		}
					// last read in the partial file section
			//		else if (readn >= (strs - 2))
			//		{
			//			end_read = reads[readn];
						// if processing last file section, the final read will end with '\0'
			//			if (file_s == file_sections - 1)
			//			{
			//				while (*end_read++ != '\0');
							// back-track 3 nucleotides (\n\0_)
			//				if (*(end_read - 2) == '\n' || *(end_read - 2) == '\r') end_read--;
							// back-track 2 nucleotides (\0_)
			//				end_read -= 2;
			//			}
						// if processing a file section > 0 and < last file section, the final read will end with '>' (beginning of split-read)
			//			else
			//			{
			//				while (*end_read++ != '>');
							// back-track 3 nucleotides (\n>_)
			//				end_read -= 3;
			//			}
			//		}
					// all other reads
			//		else end_read = reads[readn + 1] - 2;
			//	}
				// FASTQ
			//	else end_read = reads[readn] + readlen - 1;
			//	format_rev(reads[readn], end_read, &myread[0], filesig);
			////~reverse-complement read

			 // get the edit distance between reference and read
			char to_char[5] = { 'A','C','G','T','N' };
			double id = 0;
			uint32_t mismatches = 0;
			uint32_t gaps = 0;

			uint32_t ref_seq = read.hits_align_info.alignv[p].ref_seq; // ptr_alignment->ref_seq;
			//char* ref_seq_ptr = read.hits_align_info.alignv[p].ref_seq[0]; // reference_seq[(2 * ref_seq) + 1]
			//char* read_seq_ptr = myread;
			int32_t qb = read.hits_align_info.alignv[p].ref_begin1; //ptr_alignment->ref_begin1;
			int32_t pb = read.hits_align_info.alignv[p].read_begin1; //->read_begin1;
                
			Output::calcMismatchGapId(refs, read, p, mismatches, gaps, id);
			//for (uint32_t c2 = 0; c2 < ptr_alignment->cigarLen; ++c2)
			//for (uint32_t c2 = 0; c2 < read.hits_align_info.alignv[p].cigarLen; ++c2)
			//{
				//uint32_t letter = 0xf & *(ptr_alignment->cigar + c2);
				//uint32_t length = (0xfffffff0 & *(ptr_alignment->cigar + c2)) >> 4;
			//	uint32_t letter = 0xf & read.hits_align_info.alignv[p].cigar[c2];
			//	uint32_t length = (0xfffffff0 & read.hits_align_info.alignv[p].cigar[c2]) >> 4;
			//	if (letter == 0)
			//	{
			//		for (uint32_t u = 0; u < length; ++u)
			//		{
			//			if ((char)to_char[(int)*(ref_seq_ptr + qb)] != (char)to_char[(int)*(read_seq_ptr + pb)]) ++mismatches;
			//			else ++id;
			//			++qb;
			//			++pb;
			//		}
			//	}
			//	else if (letter == 1)
			//	{
			//		pb += length;
			//		gaps += length;
			//	}
			//	else
			//	{
			//		qb += length;
			//		gaps += length;
			//	}
			//}

			//int32_t align_len = abs(ptr_alignment->read_end1 + 1 - ptr_alignment->read_begin1);
			int32_t align_len = abs(read.hits_align_info.alignv[p].read_end1 + 1 - read.hits_align_info.alignv[p].read_begin1);
			int32_t total_pos = mismatches + gaps + id;
			stringstream ss;
			ss.precision(3);
			ss << (double)id / total_pos << ' ' << (double)align_len / read.hits_align_info.alignv[p].readlen;
			double align_id_round = 0.0;
			double align_cov_round = 0.0;
			ss >> align_id_round >> align_cov_round;
                
			// alignment with the highest SW score passed
			// %id and %coverage thresholds
			if ((p == index_max_score) &&
				(align_id_round >= align_id) &&
				(align_cov_round >= align_cov))
			{
				// increment number of reads passing identity
				// and coverage threshold
				readstats.total_reads_mapped_cov++;
				// do not output read for de novo OTU construction
				// (it passed the %id/coverage thresholds)
				if (de_novo_otu_gv && read.hit_denovo) read.hit_denovo != read.hit_denovo; // flip
				// fill OTU map with highest-scoring alignment for the read
				if (otumapout_gv)
				{
					// reference sequence identifier for mapped read
					//char ref_seq_arr[4000] = "";
					//char* ref_seq_arr_ptr = ref_seq_arr;
					//char* ref_seq_id_ptr = reference_seq[(2 * ref_seq)] + 1;
					//while ((*ref_seq_id_ptr != ' ') && (*ref_seq_id_ptr != '\n') && (*ref_seq_id_ptr != '\r')) 
					//	*ref_seq_arr_ptr++ = *ref_seq_id_ptr++;
					//string ref_seq_str = ref_seq_arr;
					std::string refhead = refs.buffer[read.hits_align_info.alignv[p].ref_seq].header;
					string ref_seq_str = refhead.substr(0, refhead.find(' '));

					// read identifier
					//char read_seq_arr[4000] = "";
					//char* read_seq_arr_ptr = read_seq_arr;
					//char* read_seq_id_ptr = reads[readn - 1] + 1;
					//while ((*read_seq_id_ptr != ' ') && (*read_seq_id_ptr != '\n') && (*read_seq_id_ptr != '\r')) *read_seq_arr_ptr++ = *read_seq_id_ptr++;
					//string read_seq_str = read_seq_arr;
					string read_seq_str = read.header.substr(0, read.header.find(' '));
					readstats.otu_map[ref_seq_str].push_back(read_seq_str);
				}
			}

			// output alignment to SAM or Blast-like formats
			if (samout_gv || blastout_gv)
			{
				// quality for FASTQ
				//char* read_qual = NULL;
				//if (read.format == Format::FASTA)
				//{
					// forward
					//if (ptr_alignment->strand)
				//	if (read.hits_align_info.alignv[p].strand)
				//	{
						// if second part of split read or last read in file
				//		if (((readn == 3) && (file_s > 0)) || (readn >= (strs - 2)))
				//		{
				//			read_qual = reads[readn];
				//			int8_t numnewlines = 0;
				//			while (numnewlines < 2) { if (*read_qual++ == '\n' || *read_qual++ == '\r') numnewlines++; }
				//		}
				//		else read_qual = reads[readn + 1] - read.hits_align_info.alignv[p].readlen - 1;
				//	}
					// reverse-complement
				//	else
				//	{
						// if second part of split read or last read in file
				//		if (((readn == 3) && (file_s > 0)) || (readn >= (strs - 2)))
				//		{
				//			read_qual = reads[readn];

							// last file section
				//			if (file_s == file_sections - 1)
				//			{
				//				while (*read_qual != '\0') read_qual++;
				//				if (*(read_qual - 3) == '\n' || *(read_qual - 3) == '\r') read_qual--;
								// account for '\n\0'
				//				read_qual -= 2;
				//			}
							// file section > 0 and < last file section
				//			else
				//			{
				//				while (read_qual != finalnt) read_qual++;
				//				read_qual--;
				//			}
				//		}
				//		else read_qual = reads[readn + 1] - 2;
				//	}
				//}//~if filesig == '@'

				if (blastout_gv)
				{
					uint32_t bitscore = (uint32_t)((float)((index.gumbel[index.index_num].first)
						* (read.hits_align_info.alignv[p].score1) - log(index.gumbel[index.index_num].second)) / (float)log(2));

					double evalue_score = (double)(index.gumbel[index.index_num].second) * index.full_ref[index.index_num] 
						* index.full_read[index.index_num] 
						* std::exp(-(index.gumbel[index.index_num].first) * read.hits_align_info.alignv[p].score1);
					
					Output::report_blast(output.acceptedblast, index, refs, read);
					//report_blast(
					//	output.acceptedblast, //blast output file
					//	read.hits_align_info.alignv[p], //SW alignment cigar  ptr_alignment
					//	reads[readn - 1] + 1, //read name
					//	myread, //read sequence (in integer format)
					//	read_qual, //read quality
					//	reference_seq[(2 * ref_seq)] + 1, //reference name
					//	reference_seq[(2 * ref_seq) + 1], //reference sequence
					//	evalue_score, //e-value score
					//	read.hits_align_info.alignv[p].readlen, //read length (to compute the masked regions)
					//	bitscore,
					//	ptr_alignment->strand,
					//	(double)id / total_pos, // %id
					//	(double)align_len / read.hits_align_info.alignv[p].readlen, // %query coverage
					//	mismatches,
					//	gaps
					//);
				}
				if (samout_gv)
				{
					Output::report_sam(output.acceptedblast, refs, read);
					//report_sam(
					//	acceptedsam, //sam output file
					//	ptr_alignment, //SW alignment cigar
					//	reads[readn - 1] + 1, //read name
					//	myread,//read sequence (in integer format)
					//	read_qual, //read quality
					//	reference_seq[(2 * ref_seq)] + 1, //reference name
					//	reference_seq[(2 * ref_seq) + 1], //reference sequence
					//	read.hits_align_info.alignv[p].readlen, //read length (to compute the masked regions)
					//	ptr_alignment->strand,
					//	mismatches + gaps //edit distance
					//);
				}
			}//~if (samout_gv || blastout_gv)
		}//~if alignment at current database and index part loaded in RAM
		//ptr_alignment++;
	}//~for all best alignments
} // ~writeAlignmentJob
