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
#include "kseq.h"
#include "kseq_load.hpp"
#include "traverse_bursttrie.hpp"
#include "alignment.hpp"

#include "options.hpp"
#include "ThreadPool.hpp"
#include "read.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "index.hpp"
#include "references.hpp"
#include "readsqueue.hpp"
#include "kvdb.hpp"
#include "processor.hpp"
#include "reader.hpp"
#include "writer.hpp"
#include "output.hpp"

//#ifdef _OPENMP
//#include <omp.h>
//#endif

#if defined(_WIN32)
#define O_SMR_READ_BIN O_RDONLY | O_BINARY
#else
#define O_SMR_READ_BIN O_RDONLY
#endif

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
//char complement[4] = { 3,2,1,0 };

// Callback run in a Processor thread
void parallelTraversalJob
	(
		Runopts & opts, 
		Index & index, 
		References & refs, 
		Output & output, 
		Readstats & readstats, 
		Refstats & refstats, 
		Read & read
	)
{
	read.lastIndex = index.index_num;
	read.lastPart = index.part;

	// for reverse reads
	if (!opts.forward)
	{
		// output the first num_alignments_gv alignments
		if (opts.num_alignments > 0)
		{
			// all num_alignments_gv alignments have been output
			if (read.num_alignments < 0) return;
		}
		// the maximum scoring alignment has been found, go to next read
		// (unless all alignments are being output)
		else if (opts.num_best_hits > 0 && opts.min_lis > 0 && read.max_SW_score == opts.num_best_hits)
			return;
	}

	bool read_to_count = true; // passed directly to compute_lis_alignment. TODO: What's the point?

	// find the minimum sequence length
	if (read.sequence.size() < readstats.min_read_len)
		readstats.min_read_len = static_cast<uint32_t>(read.sequence.size());

	// find the maximum sequence length
	if (read.sequence.size()  > readstats.max_read_len)
		readstats.max_read_len = static_cast<uint32_t>(read.sequence.size());

	// the read length is too short
	if (read.sequence.size()  < refstats.lnwin[index.index_num])
	{
		std::stringstream ss;
		ss << "\n  " << "\033[0;33m" << "WARNING" << COLOFF
			<< ": Processor thread: " << std::this_thread::get_id()
			<< " The read: " << read.id << " is shorter "
			<< "than " << refstats.lnwin[index.index_num] << " nucleotides, by default it will not be searched\n";
		std::cout << ss.str(); ss.str("");

		read.isValid = false;
		return;
	}

	uint32_t windowshift = opts.skiplengths[index.index_num][0];
	// keep track of windows which have been already traversed in the burst trie
	vector<bool> read_index_hits(read.sequence.size());

	uint32_t pass_n = 0; // Pass number (possible value 0,1,2)
	uint32_t max_SW_score = read.sequence.size() *opts.match; // the maximum SW score attainable for this read

	std::vector<MYBITSET> vbitwindowsf;
	std::vector<MYBITSET> vbitwindowsr;

	// TODO: these 2 lines need to be calculated once per index part. Move to index or some new calculation object
	uint32_t bit_vector_size = (refstats.partialwin[index.index_num] - 2) << 2; // on each index part. Init in Main, accessed in Worker
	uint32_t offset = (refstats.partialwin[index.index_num] - 3) << 2; // on each index part

	uint32_t minoccur = 0; // TODO: never updated. Always 0. What's the point?

	// loop for each new Pass to granulate seed search intervals
	for (bool search = true; search; )
	{
		uint32_t numwin = (read.sequence.size()
			- refstats.lnwin[index.index_num]
			+ windowshift) / windowshift; // number of k-mer windows fit along the sequence

		uint32_t win_index = 0; // index of the window's first char in the sequence e.g. 0, 18, 36 if window.length = 18
		// iterate the windows
		for (uint32_t win_num = 0; win_num < numwin; win_num++)
		{
			// skip position, seed at this position has already been searched for in a previous Pass
			//if (read_index_hits[win_index]) goto check_score;
			// search position, set search bit to true
			//else read_index_hits[win_index].flip();
			if (!read_index_hits[win_index])
			{
				read_index_hits[win_index].flip();
				// this flag it set to true if a match is found during
				// subsearch 1(a), to skip subsearch 1(b)
				bool accept_zero_kmer = false;
				// ids for k-mers that hit the database
				vector<id_win> id_hits; // TODO: why not to add directly to 'id_win_hits'? - because id_win_hits may contain hits from different index parts.
				vbitwindowsf.resize(bit_vector_size);
				std::fill(vbitwindowsf.begin(), vbitwindowsf.end(), 0);

				init_win_f(&read.isequence[win_index + refstats.partialwin[index.index_num]],
					&vbitwindowsf[0],
					&vbitwindowsf[4],
					refstats.numbvs[index.index_num]);

				uint32_t keyf = 0;
				char *keyf_ptr = &read.isequence[win_index];
				// build hash for first half windows (foward and reverse)
				// hash is just a numeric value formed by the chars of a string consisting of '0','1','2','3'
				// e.g. "2233012" -> b10.1011.1100.0110 = x2BC6 = 11206
				for (uint32_t g = 0; g < refstats.partialwin[index.index_num]; g++)
				{
					(keyf <<= 2) |= (uint32_t)(*keyf_ptr);
					++keyf_ptr;
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
						refstats.partialwin[index.index_num],
						opts
					);
				} //~if exact half window exists in the burst trie

				// only search if an exact match has not been found
				if (!accept_zero_kmer)
				{
					vbitwindowsr.resize(bit_vector_size);
					std::fill(vbitwindowsr.begin(), vbitwindowsr.end(), 0);

					// build the first bitvector window
					init_win_r(&read.isequence[win_index + refstats.partialwin[index.index_num] - 1],
						&vbitwindowsr[0],
						&vbitwindowsr[4],
						refstats.numbvs[index.index_num]);

					uint32_t keyr = 0;
					char *keyr_ptr = &read.isequence[win_index + refstats.partialwin[index.index_num]];

					// build hash for first half windows (foward and reverse)
					for (uint32_t g = 0; g < refstats.partialwin[index.index_num]; g++)
					{
						(keyr <<= 2) |= (uint32_t)(*keyr_ptr); //  - '0'
						++keyr_ptr;
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
							&vbitwindowsr[0], /* win1r_k1_ptr */
							&vbitwindowsr[offset], /* win1r_k1_full */
							accept_zero_kmer,
							id_hits,
							read.id,
							win_index,
							refstats.partialwin[index.index_num], 
							opts);
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

			// continue read analysis if threshold seeds were matched
			if (win_num == numwin - 1)
			{
				compute_lis_alignment(
					read, opts, index, refs, readstats, refstats, output,
					search, // returns False if the alignment is found -> stop searching
					max_SW_score,
					read_to_count
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
							(opts.skiplengths[index.index_num][pass_n] == opts.skiplengths[index.index_num][pass_n + 1])) ++pass_n;
						if (++pass_n > 2) search = false;
						// set interval skip length for next Pass
						else windowshift = opts.skiplengths[index.index_num][pass_n];
					}
				}
				break; // do not offset final window on read
			}//~( win_num == NUMWIN-1 )
			win_index += windowshift;
		}//~for (each window)                
			//~while all three window skip lengths have not been tested, or a match has not been found
	}// ~while (search);

	// the read didn't align (for --num_alignments [INT] option),
	// output null alignment string
	if (!read.hit && !opts.forward && (opts.num_alignments > -1))
	{
		// do not output read for de novo OTU clustering
		// (it did not pass the E-value threshold)
		if (opts.de_novo_otu && read.hit_denovo) read.hit_denovo = !read.hit_denovo; // flip
	}//~if read didn't align

	if (opts.de_novo_otu && read.hit_denovo)
		++readstats.total_reads_denovo_clustering;
} // ~parallelTraversalJob


void paralleltraversal(Runopts & opts)
{
	std::stringstream ss;

	unsigned int numCores = std::thread::hardware_concurrency(); // find number of CPU cores
	ss << "CPU cores on this machine: " << numCores << std::endl; // 8
	std::cout << ss.str(); ss.str("");

	// Init thread pool with the given number of threads
	int numThreads = 2 * opts.num_fread_threads + opts.num_proc_threads;
	if (numThreads > numCores) {
		ss << "WARN: Number of cores: " << numCores << " is less than number allocated threads " << numThreads << std::endl;
		std::cout << ss.str(); ss.str("");
	}

	ThreadPool tpool(numThreads);
	KeyValueDatabase kvdb(opts.kvdbPath);
	ReadsQueue readQueue("read_queue", QUEUE_SIZE_MAX, 1); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", QUEUE_SIZE_MAX, opts.num_proc_threads); // shared: Processor pushes, Writer pops
	Readstats readstats(opts);
	Refstats refstats(opts, readstats);
	Output output(opts, readstats);
	Index index;
	References refs;

	int loopCount = 0; // counter of total number of processing iterations

	// perform alignment
	auto starts = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;
	// loop through every index passed to option --ref (ex. SSU 16S and SSU 18S)
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); ++index_num)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
		{
			ss << "\tLoading index " << index_num << " part " << idx_part + 1 << "/" << refstats.num_index_parts[index_num] << " ... ";
			std::cout << ss.str(); ss.str("");

			starts = std::chrono::high_resolution_clock::now();
			index.load(index_num, idx_part, opts, refstats);
			refs.load(index_num, idx_part, opts, refstats);
			elapsed = std::chrono::high_resolution_clock::now() - starts; // ~20 sec Debug/Win
//			std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - t);

			ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << "] sec\n";
			std::cout << ss.str(); ss.str("");

			starts = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < opts.num_fread_threads; i++)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), opts, readQueue, kvdb, loopCount));
				tpool.addJob(Writer("writer_" + std::to_string(i), writeQueue, kvdb));
			}

			// add processor jobs
			for (int i = 0; i < opts.num_proc_threads; i++)
			{
				tpool.addJob(Processor("proc_" + std::to_string(i), readQueue, writeQueue, opts, index, refs, output, readstats, refstats, parallelTraversalJob));
			}
			++loopCount;

			tpool.waitAll(); // wait till all reads are processed against the current part
			index.clear();
			refs.clear();
			writeQueue.reset(opts.num_proc_threads);
			readQueue.reset(1);

			elapsed = std::chrono::high_resolution_clock::now() - starts;
			ss << "    Done index " << index_num << " Part: " << idx_part + 1 
				<< " Time: " << std::setprecision(2) << std::fixed << elapsed.count() << " sec\n";
			std::cout << ss.str(); ss.str("");
		} // ~for(idx_part)
	} // ~for(index_num)
	// store readstats calculated in alignment
	kvdb.put("Readstats", readstats.toString());
} // ~paralleltraversal
