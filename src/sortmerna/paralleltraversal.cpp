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


#if defined(_WIN32)
#define O_SMR_READ_BIN O_RDONLY | O_BINARY
#else
#define O_SMR_READ_BIN O_RDONLY
#endif

// forward
int clear_dir(std::string dpath);

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

/* 
 * Callback run in a Processor thread
 * Called on each index * index_part * read.num_strands
 *
 * @param isLastStrand Boolean flags when the last strand is passed for matching
 */
void alignmentCb
	(
		Runopts & opts, 
		Index & index, 
		References & refs, 
		Output & output, 
		Readstats & readstats, 
		Refstats & refstats, 
		Read & read,
		bool isLastStrand
	)
{
	read.lastIndex = index.index_num;
	read.lastPart = index.part;

	// for reverse reads
	if (read.reversed)
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
		ss << "\n  " << YELLOW << "WARNING" << COLOFF
			<< __FILE__ << ":" << __LINE__ << ": Processor thread: " << std::this_thread::get_id()
			<< " The read: " << read.id << " is shorter than "
			<< refstats.lnwin[index.index_num] << " nucleotides, by default it will not be searched\n";
		std::cout << ss.str(); ss.str("");

		read.isValid = false;
		return;
	}

	uint32_t windowshift = opts.skiplengths[index.index_num][0];
	// keep track of windows (read positions) which have been already traversed in the burst trie
	// initially all False
	vector<bool> read_pos_searched(read.sequence.size());

	uint32_t pass_n = 0; // Pass number (possible value 0,1,2)
	uint32_t max_SW_score = read.sequence.size() *opts.match; // the maximum SW score attainable for this read

	std::vector<UCHAR> bitvec; // window (prefix/suffix) bitvector

	// TODO: below 2 values are unique per index part. Move to index?
	uint32_t bitvec_size = (refstats.partialwin[index.index_num] - 2) << 2; // e.g. 9 - 2 = 0000 0111 << 2 = 0001 1100 = 28
	// Does this mark where in 32-bit the bitvector starts?
	uint32_t offset = (refstats.partialwin[index.index_num] - 3) << 2; // e.g. 9 - 3 = 0000 0110 << 2 = 0001 1000 = 24

	// loop search positions on the read in multiple passes
	// changing the step (windowshift) when necessary
	for (bool search = true; search; )
	{
		// number of k-mer windows fit along the read given 
		// the window size and a search step (windowshift)
		uint32_t numwin = ( read.sequence.size()
				- refstats.lnwin[index.index_num]
				+ windowshift ) / windowshift;

		uint32_t win_pos = 0; // position (index) of the window's first char on the read i.e. [0..read.sequence.length-1]
		// iterate the windows
		for (uint32_t win_num = 0; win_num < numwin; win_num++)
		{
			if (read.is04) read.flip34(); // Make sure the read is in 03 encoding for index search

			// skip position when the seed at this position has already been searched for in a previous Passes
			if (!read_pos_searched[win_pos])
			{
				read_pos_searched[win_pos].flip(); // mark position as searched
				// this flag it set to true if a match is found during
				// subsearch 1(a), to skip subsearch 1(b)
				bool accept_zero_kmer = false;
				// ids for k-mers that hit the database
				vector<id_win> id_hits; // TODO: add directly to 'id_win_hits'? - No, id_win_hits may contain hits from different index parts.

				bitvec.resize(bitvec_size);
				std::fill(bitvec.begin(), bitvec.end(), 0);

				init_win_f(&read.isequence[win_pos + refstats.partialwin[index.index_num]],
					&bitvec[0],
					&bitvec[4],
					refstats.numbvs[index.index_num]);

				// the hash of the first half of the kmer window
				uint32_t keyf = read.hashKmer(win_pos, refstats.partialwin[index.index_num]);

				// TODO: remove in production
				if (index.lookup_tbl.size() <= keyf) {
					std::stringstream ss;
					size_t vsize = index.lookup_tbl.size();
					uint16_t idxn = index.index_num;
					uint16_t idxp = index.part;
					unsigned int id = read.id;
					bool is03 = read.is03;
					bool is04 = read.is04;
					ss << __FILE__ << ":" << __LINE__
						<< " ERROR: lookup index: " << keyf << " is larger than lookup_tbl.size: " << vsize 
						<< " Index: " << idxn
						<< " Part: " << idxp
						<< " Read.id: " << id
						<< " Read.is03: " << is03
						<< " Read.is04: " << is04
						<< " Aborting.." << std::endl;
					std::cout << ss.str();
					exit(EXIT_FAILURE);
				}

				// do traversal if the exact half window exists in the burst trie
				if ( index.lookup_tbl[keyf].count > opts.minoccur && index.lookup_tbl[keyf].trie_F != NULL )
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
						&bitvec[0],
						&bitvec[offset],
						accept_zero_kmer,
						id_hits,
						read.id,
						win_pos,
						refstats.partialwin[index.index_num],
						opts
					);
				} //~if exact half window exists in the burst trie

				// only search reversed kmer if an exact match has not been found for the forward
				if (!accept_zero_kmer)
				{
					//bitvec.resize(bitvec_size);
					std::fill(bitvec.begin(), bitvec.end(), 0);

					// init the first bitvector window
					init_win_r(&read.isequence[win_pos + refstats.partialwin[index.index_num] - 1],
						&bitvec[0],
						&bitvec[4],
						refstats.numbvs[index.index_num]);

					// the hash of the second (rear) half of the kmer window
					uint32_t keyr = read.hashKmer(win_pos + refstats.partialwin[index.index_num], refstats.partialwin[index.index_num]);

					// TODO: remove in production
					if (index.lookup_tbl.size() <= keyr) {
						std::stringstream ss;
						size_t vsize = index.lookup_tbl.size();
						uint16_t idxn = index.index_num;
						uint16_t idxp = index.part;
						unsigned int id = read.id;
						bool is03 = read.is03;
						bool is04 = read.is04;
						ss << __LINE__ << " Thread: " << std::this_thread::get_id()
							<< " ERROR: lookup index: " << keyr << " is larger than lookup_tbl.size: " << vsize
							<< " Index: " << idxn
							<< " Part: " << idxp
							<< " Read.id: " << id
							<< " Read.is03: " << is03
							<< " Read.is04: " << is04
							<< " Aborting.." << std::endl;
						std::cout << ss.str();
						exit(EXIT_FAILURE);
					}

					// continue subsearch (1)(b)
					if ( index.lookup_tbl[keyr].count > opts.minoccur && index.lookup_tbl[keyr].trie_R != NULL )
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
							&bitvec[0],
							&bitvec[offset],
							accept_zero_kmer,
							id_hits,
							read.id,
							win_pos,
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
			} // ~if not read_pos_searched[win_pos]

			// continue read analysis if threshold seeds were matched
			if (win_num == numwin - 1)
			{
				compute_lis_alignment(
					read, opts, index, refs, readstats, refstats,
					search, // returns False if the alignment is found -> stop searching
					max_SW_score,
					read_to_count
				);

				// the read was not accepted at current window shift,
				// use the next (smaller) window shift
				if (search)
				{
					// last (3rd) Pass has been made
					if (pass_n == 2) search = false;
					else
					{
						// the next interval size equals to the current one, skip it
						while ( pass_n < 3 &&
							opts.skiplengths[index.index_num][pass_n] == opts.skiplengths[index.index_num][pass_n + 1] )
							++pass_n;
						if (++pass_n > 2) search = false;
						// set interval skip length for next Pass
						else windowshift = opts.skiplengths[index.index_num][pass_n];
					}
				}
				break; // last possible position reached for given window and skip length -> go to the next skip length
			}//~( win_num == NUMWIN-1 )
			win_pos += windowshift;
		}//~for (each window)                
			//~while all three window skip lengths have not been tested, or a match has not been found
	}// ~while (search);

	// the read didn't align (for --num_alignments [INT] option),
	// output null alignment string
	if (isLastStrand && !read.hit && opts.num_alignments > -1) // !opts.forward
	{
		// do not output read for de novo OTU clustering
		// (it did not pass the E-value threshold)
		if (opts.de_novo_otu) read.hit_denovo = false;
	}//~if read didn't align
} // ~alignmentCb

// called from main
void align(Runopts & opts, Readstats & readstats, Output & output)
{
	std::stringstream ss;

	unsigned int numCores = std::thread::hardware_concurrency(); // find number of CPU cores

	// Init thread pool with the given number of threads
	int numProcThread = 0;
	if (opts.num_proc_thread == 0) {
		numProcThread = numCores; // default
		ss << "paralleltraversal: Using default number of Processor threads = num CPU cores: " << numCores << std::endl; // 8
		std::cout << ss.str(); ss.str("");
	}
	else
	{
		numProcThread = opts.num_proc_thread; // set using '--thread'
		ss << "paralleltraversal: Using number of Processor threads set in run options: " << numProcThread << std::endl; // 8
		std::cout << ss.str(); ss.str("");
	}

	int numThreads = opts.num_read_thread + opts.num_write_thread + numProcThread;

	ss << "Number of cores: " << numCores 
		<< " Read threads:  " << opts.num_read_thread
		<< " Write threads: " << opts.num_write_thread
		<< " Processor threads: " << numProcThread
		<< std::endl;
	std::cout << ss.str(); ss.str("");

	ThreadPool tpool(numThreads);
	clear_dir(opts.kvdbPath);
	KeyValueDatabase kvdb(opts.kvdbPath);
	ReadsQueue readQueue("read_queue", QUEUE_SIZE_MAX, opts.num_read_thread); // shared: Processor pops, Reader pushes
	ReadsQueue writeQueue("write_queue", QUEUE_SIZE_MAX, numProcThread); // shared: Processor pushes, Writer pops
	Refstats refstats(opts, readstats);
	Index index;
	References refs;

	int loopCount = 0; // counter of total number of processing iterations

	// perform alignment
	auto starts = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;

	// loop through every index passed to option '--ref'
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); ++index_num)
	{
		// iterate every part of an index
		for (uint16_t idx_part = 0; idx_part < refstats.num_index_parts[index_num]; ++idx_part)
		{
			ss << std::endl << "    Loading index " << index_num << " part " << idx_part + 1 << "/" << refstats.num_index_parts[index_num] << " ... ";
			std::cout << ss.str(); ss.str("");

			starts = std::chrono::high_resolution_clock::now();
			index.load(index_num, idx_part, opts, refstats);
			refs.load(index_num, idx_part, opts, refstats);
			elapsed = std::chrono::high_resolution_clock::now() - starts; // ~20 sec Debug/Win
//			std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - t);

			ss << "done [" << std::setprecision(2) << std::fixed << elapsed.count() << "] sec\n";
			std::cout << ss.str(); ss.str("");

			starts = std::chrono::high_resolution_clock::now();
			for (int i = 0; i < opts.num_read_thread; i++)
			{
				tpool.addJob(Reader("reader_" + std::to_string(i), opts, readQueue, kvdb, loopCount));
			}

			for (int i = 0; i < opts.num_write_thread; i++)
			{
				tpool.addJob(Writer("writer_" + std::to_string(i), writeQueue, kvdb));
			}

			// add processor jobs
			for (int i = 0; i < numProcThread; i++)
			{
				tpool.addJob(Processor("proc_" + std::to_string(i), readQueue, writeQueue, opts, index, refs, output, readstats, refstats, alignmentCb));
			}
			++loopCount;

			tpool.waitAll(); // wait till all reads are processed against the current part
			index.clear();
			refs.clear();
			writeQueue.reset(numProcThread);
			readQueue.reset(opts.num_read_thread);

			elapsed = std::chrono::high_resolution_clock::now() - starts;
			ss << "    paralleltraversal: Done index " << index_num << " Part: " << idx_part + 1 
				<< " Time: " << std::setprecision(2) << std::fixed << elapsed.count() << " sec" << std::endl;
			std::cout << ss.str(); ss.str("");
		} // ~for(idx_part)
	} // ~for(index_num)

	// store readstats calculated in alignment
	kvdb.put("Readstats", readstats.toString());
} // ~align
