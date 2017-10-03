/**
 * @file paralleltraversal.hpp
 * @brief Function and variable definitions for paralleltraversal.cpp
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
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 */

#pragma once

#include <iomanip>
#include <string>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <chrono>
 
#include "outputformats.hpp"
#include "load_index.hpp"
#include "traverse_bursttrie.hpp"
#include "alignment.hpp"
#include "mmap.hpp"
#include "kseq_load.hpp"
#include "options.hpp"


 /*! @fn check_file_format()
	 @brief check input reads file format (FASTA, FASTQ or unrecognized)
	 @param char* inputreads
	 @param char& filesig
	 @return bool
	 @version Feb 15, 2016
  */
bool
check_file_format(
	char* inputreads /**< pointer to query reads file */,
	char& filesig /**< first character of sequence label */);

/*! @fn compute_read_stats()
	@brief compute total number of reads in file and their combined length
	@param char* inputreads
	@param uint64_t& number_total_read
	@param uint64_t& full_read_main
	@param off_t& full_file_size
	@version Feb 15, 2016
 */
void compute_read_stats(
	char* inputreads /**< pointer to query reads file */,
	uint64_t& number_total_read /**< total number of reads */,
	uint64_t& full_read_main /**< total number of nucleotides in all reads */,
	off_t& full_file_size /**< the size of the full reads file (in bytes) */);


/*! @fn paralleltraversal()
	@brief Traverse the query input and indexed database and output
		   alignments passing the E-value threshold
	@detail The following methods will be executed:
	<ol>
	  <li> divide large read files into mmap'd regions,
		   taking into account the read (and its pair) which may
		   be split between two file sections </li>
	  <li> compute the gumbel parameters (lamda and K) using ALP,
		   load the index fully or in parts (depending on how
		   it was built) </li>
	  <li> using 3 intervals, scan over the read and collect all
		   L-mers on the read which match to the reference
		   index with at most 1 error. This is done using
		   parallel traversal between the index and the
		   Levenshtein automaton </li>
	  <li> if enough L-mers were collected, extend them into
		   longer matches using the Longest Increasing
		   subsequence (LIS) of positions where the L-mers
		   matched on the read and the associated reference
		   sequences </li>
	  <li> if the LIS is long enough, use the starting positions
		   of the LIS to estimate the starting position
		   of an alignment and pass this reference segment and
		   read to SSW </li>
	  <li> if the alignment score is at least the minimum score
		   corresponding to the E-value threshold, keep the read,
		   otherwise continue searching for other LIS or more
		   L-mers using smaller intervals </li>
	</ol>

	@param char* inputreads
	@param bool have_reads_gz
	@param *ptr_filetype_ar
	@param *ptr_filetype_or
	@param long match
	@param long mismatch
	@param long gap_open
	@param long gap_extension
	@param long score_N
	@param vector< vector<uint32_t> >
	@param int argc
	@param char **argv
	@param bool yes_SQ
	@param vector< pair<string,string> >& myfiles
	@return void
	@version Feb 10, 2016
*/
void paralleltraversal(
	char* inputreads /**< pointer to query reads file */,
	bool have_reads_gz /**< if true, input reads file is in compressed format */,
	char* ptr_filetype_ar /**< pointer to string for aligned seqeunces filepath */,
	char* ptr_filetype_or /**< pointer to string for rejected sequences filepath */,
	long match /**< SW match reward score (positive) */,
	long mismatch /**< SW mismatch penalty score (negative) */,
	long gap_open /**< SW gap open penalty score (positive) */,
	long gap_extension /**< SW gap extend penalty score (positive) */,
	long score_N /**< SW penalty for ambiguous nucleotide (negative) */,
	vector< vector<uint32_t> >& skiplengths /**< skiplengths, three intervals at which to place seeds on read */,
	int argc /**< number of arguments passed to SortMeRNA */,
	char **argv /**< argument string passed to SortMeRNA */,
	bool yes_SQ /**< if true, include @SQ tags in SAM output */,
	vector< pair<string, string> >& myfiles /**< vector of (FASTA file, index name) pairs for loading index */,
	bool exit_early /**< if true, exit program if reads file is not FASTA or FASTQ, or reads files or reference file is empty */);

void paralleltraversal2(Runopts & opts);

//
// "Producer-Consumer" concurrent pattern participants
//

const char FASTA_HEADER_START = '>';
const char FASTQ_HEADER_START = '@';
const char WIN_ENDLINE = '\r';
const char LIN_ENDLINE = '\n';
const int QUEUE_SIZE_MAX = 10; // TODO: set through process options
const int NUM_PROC_THREADS = 3; // Default number of reads processor threads. Change through process options.

// Collective Statistics for all Reads. Encapsulates old 'compute_read_stats' logic and results
struct ReadStatsAll {
	// Synchronized - Compute in worker thread (per read) - shared
	uint32_t min_read_len = READLEN; // length of the shortest Read in the Reads file. 
	uint32_t max_read_len; // length of the longest Read in the Reads file.
	uint64_t total_reads_mapped; // total number of reads mapped passing E-value threshold. Computed in 'compute_lis_alignment' in a worker thread i.e. per read.

	// TODO: move to Readrec and get rid of this vector
	//std::vector<bool> read_hits; // flags if a read was aligned i.e. match found. Each value represents a read. True when the read was matched/aligned.

	// TODO: move to Readrec and get rid of this vector
	// bits representing all reads. An accepted read with < %id and < %coverage is set to false (0)
	//std::vector<bool> read_hits_denovo;

	// TODO: move to Readrec and get rid of this vector
	// array of uint16_t to represent all reads, if the read was aligned with a maximum SW score, its number of alignments is incremeted by 1
	//uint16_t *read_max_SW_score; // SW (Smith-Waterman) Max SW score of a read -> Readrec.max_SW_score
	
	// Non-synchronized - Compute once by calculate
	uint64_t number_total_read; // total number of reads in file.
	off_t    full_file_size; // the size of the full reads file (in bytes).
	uint64_t full_read_main; // total number of nucleotides in all reads.

	ReadStatsAll() {}
	ReadStatsAll(std::string & readsfile) {}
	
	~ReadStatsAll() {}

	// calculate statistics from readsfile see "compute_read_stats"
	void calculate(std::string & readsfile) {}
	// called from Main thread once when 'number_total_read' is known
	//void set_read_hits() {
	//	read_hits.resize(2*number_total_read); // why twice the number of reads? Because original **reads array has 2 lines per read: header and sequence.
	//	std::fill(read_hits.begin(), read_hits.end(), false);
	//}
	//void set_read_hits_denovo() {
	//	read_hits_denovo.resize(2*number_total_read);
	//	std::fill(read_hits_denovo.begin(), read_hits_denovo.end(), true);
	//}
	//void set_read_max_SW_score() {
	//	read_max_SW_score = new uint16_t[2*number_total_read];
	//	memset(read_max_SW_score, 0, sizeof(uint16_t)*(2*number_total_read));
	//}
};

// Match (search) results for a particular Read
struct MatchResults {
	bool hit = false; // indicates that a match for this Read has been found
	bool hit_denovo = true;
	bool null_align_output = false; // flags NULL alignment was output to file (needs to be done once only)
	uint16_t max_SW_score; // Max Smith-Waterman score
	int32_t num_alignments; // number of alignments to output per read
	alignment_struct hits_align_info; // 

	MatchResults() {}
	~MatchResults() {}
};

/**
 * Wrapper of a Reads' record and its Match results
 */
struct Read {
	std::string header;
	std::string sequence;
	std::string quality; // "" (fasta) | "xxx..." (fastq)
	std::string format = "fasta"; // fasta | fastq

	// calculated
	std::string sequenceInt; // sequence in Integer alphabet: [A,C,G,T] -> [0,1,2,3]
	std::vector<int> ambiguous_nt; // positions of ambiguous nucleotides in the sequence (as defined in nt_table/load_index.cpp)

	MatchResults matches;

	Read() {}
	Read(std::string header, std::string sequence, std::string quality, std::string format) :
		header(header), sequence(sequence), quality(quality), format(format) 
	{
		validate();
		sequenceToInt();
	}

	~Read() {}

	// convert sequence to "sequenceInt" and populate "ambiguous_nt"
	void sequenceToInt() {
		for (std::string::iterator it = sequence.begin(); it != sequence.end(); ++it)
		{
			char c = (4 == nt_table[(int)*it]) ? 0 : nt_table[(int)*it];
			sequenceInt.append(1, nt_table[(int)*it]);
			if (c == 0) { // ambiguous nt
				ambiguous_nt.push_back(sequenceInt.size() - 1); // i.e. add current position to the vector
			}
		}
	}

	void validate() {
		if (sequence.size() > READLEN)
		{
			fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] at least one of your reads is > %d nt \n",
				startColor, "\033[0m", __LINE__, __FILE__, READLEN);
			fprintf(stderr, "  Please check your reads or contact the authors.\n");
			exit(EXIT_FAILURE);
		}
	} // ~validate
}; // ~struct Read

/**
 * Qeueu for Reads' records. Concurrently accessed by the Reader (producer) and the Processors (consumers)
 */
struct ReadsQueue {
	std::queue<Read> recs;
	int capacity; // max size of the queue
	int queueSizeAvr; // average size of the queue
	bool doneAdding; // flag indicating no more records will be added
	std::mutex lock;
	std::condition_variable cv;
	int numPushed = 0;
	int numPopped = 0;

	ReadsQueue(int capacity) : capacity(capacity), doneAdding(false) {}

	~ReadsQueue() {}

	void push(Read & readsrec) {
		std::unique_lock<std::mutex> l(lock);
		cv.wait(l, [this] {return recs.size() < capacity;});
		recs.push(readsrec);
//		l.unlock();
		++numPushed;
		cv.notify_one();
	}

	Read pop() {
		std::unique_lock<std::mutex> l(lock);
		cv.wait(l, [this] { return doneAdding || !recs.empty();}); //  if False - keep waiting, else - proceed.
		Read rec;
		if (!recs.empty()) {
			rec = recs.front();
			recs.pop();
		}
//		l.unlock(); // probably redundant. The lock will be released when function returns
		++numPopped;
		if (numPopped % 10000 == 0)
			printf("\rPushed: %d Popped: %d", numPushed, numPopped);
		cv.notify_one();
		return rec;
	}

	void mDoneAdding() {
		std::lock_guard<std::mutex> l(lock);
		doneAdding = true;
		cv.notify_one(); // otherwise pop can stuck not knowing the adding stopped
	}
}; // ~struct ReadsQueue

// Queue for ReadStats objects ready to be written to disk
struct WriteQueue {
	WriteQueue(){}
	~WriteQueue() {}
};

// Queue for Reads IDs used to synchronize ordering of the records in Reads file and ReadStats file
struct ReadWriterCounterQueue {
	ReadWriterCounterQueue(){}
	~ReadWriterCounterQueue() {}
};

// reads Reads and ReadStats files, generates Read objects and pushes them onto ReadsQueue
class Reader {
public:
	Reader(int id, ReadsQueue & readsQueue, std::string & readsfile) : id(id), readsQueue(readsQueue), readsfile(readsfile) {}
	void operator()() { read(); }
	void read();
private:
	int id;
	ReadsQueue & readsQueue;
	std::string & readsfile;
};

class Processor {
public:
	Processor(int id, 
		ReadsQueue & recs,
		ReadStatsAll & readstats,
		Index & index, 
		std::function<void(Read)> callback
	) :
		id(id), 
		recs(recs),
		readstats(readstats),
		index(index),
		callback(callback) {}

	void operator()() { process(); }
	void process() {
		Read rec;
		int count = 0;
//		std::cout << "Processor " << id << " (" << std::this_thread::get_id() << ") started" << std::endl;
		printf("Processor %d (%d) started\n", id, std::this_thread::get_id());
		for (;;) {
			rec = recs.pop();
			// std::cout << "Processor " << id << " Popped: " << rec.header << std::endl;
			callback(rec);
			if (rec.header == "") break;
			count++;
		}
//		std::cout << "Processor " << id << " (" << std::this_thread::get_id() << ") done. Processed " << count << std::endl;
		printf("Processor %d (%d) done. Processed %d reads\n", id, std::this_thread::get_id());
	}
private:
	int id;
	ReadsQueue & recs;
	ReadStatsAll & readstats;
	Index & index;
	std::function<void(Read)> callback;
};

class Writer {
	Writer(){}
	~Writer() {}
};

// ~PARALLELTRAVERSAL_H
