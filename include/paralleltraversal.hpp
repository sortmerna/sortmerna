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

#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

 /*! @fn check_file_format()
	 @brief check input reads file format (FASTA, FASTQ or unrecognized)
	 @param char* inputreads
	 @param char& filesig
	 @return bool
	 @version Feb 15, 2016
  */
bool check_file_format(
	const char* inputreads /**< pointer to query reads file */,
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
struct ReadStats {
	Runopts & opts;

	// Synchronized - Compute in worker thread (per read) - shared
	uint32_t min_read_len = READLEN; // length of the shortest Read in the Reads file. 
	uint32_t max_read_len; // length of the longest Read in the Reads file.
	uint64_t total_reads_mapped; // total number of reads mapped passing E-value threshold. Computed in 'compute_lis_alignment' in a worker thread i.e. per read.
	char filesig = '>';

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

	ReadStats(Runopts & opts): opts(opts) {}
	
	~ReadStats() {}

	// calculate statistics from readsfile see "compute_read_stats"
	void calculate();
	bool check_file_format();
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


/**
 * Wrapper of a Reads' record and its Match results
 */
struct Read {
	int id = 0; // number of the read in the reads file

	std::string header;
	std::string sequence;
	std::string quality; // "" (fasta) | "xxx..." (fastq)
	std::string format = "fasta"; // fasta | fastq

	// calculated
	std::string sequenceInt; // sequence in Integer alphabet: [A,C,G,T] -> [0,1,2,3]
	std::vector<int> ambiguous_nt; // positions of ambiguous nucleotides in the sequence (as defined in nt_table/load_index.cpp)

	// matching results
	bool hit = false; // indicates that a match for this Read has been found
	bool hit_denovo = true;
	bool null_align_output = false; // flags NULL alignment was output to file (needs to be done once only)
	uint16_t max_SW_score = 0; // Max Smith-Waterman score
	int32_t num_alignments = 0; // number of alignments to output per read
	alignment_struct hits_align_info;

	Read(): hits_align_info(0,0,0,0,0) {
		// create new instance of alignments
		//hits_align_info.max_size = 0;
		//hits_align_info.size = 0;
		//hits_align_info.min_index = 0;
		//hits_align_info.max_index = 0;
		//hits_align_info.ptr = new s_align[1](); // see alignment.cpp
		//hits_align_info.ptr->cigar = 0;
		//hits_align_info.ptr->cigar = new uint32_t[1];
		//hits_align_info.ptr->cigarLen = 0;
		//hits_align_info.ptr->index_num = 0;
		//hits_align_info.ptr->part = 0;
		//hits_align_info.ptr->readlen = 0;
		//hits_align_info.ptr->read_begin1 = 0;
		//hits_align_info.ptr->read_end1 = 0;
		//hits_align_info.ptr->ref_begin1 = 0;
		//hits_align_info.ptr->ref_end1 = 0;
		//hits_align_info.ptr->ref_seq = 0;
		//hits_align_info.ptr->score1 = 0;
		//hits_align_info.ptr->strand = 0;
	}

	Read(int id, std::string header, std::string sequence, std::string quality, std::string format) :
		id(id), header(header), sequence(sequence), quality(quality), format(format) 
	{
		validate();
		sequenceToInt();
	}

	~Read() { 
		if (hits_align_info.ptr != 0) {
			//	delete[] read.hits_align_info.ptr->cigar;
			//	delete read.hits_align_info.ptr;
			delete hits_align_info.ptr; 
		} 
	}

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


	std::string matchesToJson() {
		rapidjson::StringBuffer sbuf;
		rapidjson::Writer<rapidjson::StringBuffer> writer(sbuf);

		writer.StartObject();
		writer.Key("hit");
		writer.Bool(hit);
		writer.Key("hit_denovo");
		writer.Bool(hit_denovo);
		writer.Key("null_align_output");
		writer.Bool(null_align_output);
		writer.Key("max_SW_score");
		writer.Uint(max_SW_score);
		writer.Key("num_alignments");
		writer.Int(num_alignments);

		writer.Key("hits_align_info");
		writer.StartObject();
		writer.String("max_size");
		writer.Uint(10);
		writer.EndObject();

		writer.EndObject();

		return sbuf.GetString();
	} // ~Read::matchesToJsonString

	std::string matchesToString() {
		if (hits_align_info.ptr == 0)
			return "";

		// total size of this object
		const uint32_t size =
			3 * sizeof(bool)
			+ sizeof(uint16_t)
			+ sizeof(int32_t)
			/* alignment_struct */
			+ 4 * sizeof(uint32_t)
			/* s_align */
			+ hits_align_info.size * (
				hits_align_info.ptr->cigarLen * (sizeof(uint32_t))
				+ sizeof(uint32_t)
				+ 4 * sizeof(int32_t)
				+ sizeof(uint32_t)
				+ 4 * sizeof(uint16_t)
				+ sizeof(bool));

		char *target = new char[size];
		uint32_t offset = 0;
		memcpy(target, (void*)&hit, 3 * sizeof(bool));
		offset += 3 * sizeof(bool);
		memcpy(target + offset, (void*)&max_SW_score, sizeof(uint16_t));
		offset += sizeof(uint16_t);
		memcpy(target + offset, (void*)&num_alignments, sizeof(int32_t));
		offset += sizeof(int32_t);
		memcpy(target + offset, (void*)&hits_align_info, 4 * sizeof(uint32_t));
		offset += 4 * sizeof(uint32_t);
		// cigar
		memcpy(target + offset, (void*)&hits_align_info.ptr->cigar, hits_align_info.ptr->cigarLen * (sizeof(uint32_t)));
		offset += hits_align_info.ptr->cigarLen * (sizeof(uint32_t));
		memcpy(target + offset, (void*)&hits_align_info.ptr->ref_seq, sizeof(uint32_t));
		offset += sizeof(uint32_t);
		memcpy(target + offset, (void*)&hits_align_info.ptr->ref_begin1, 4 * sizeof(int32_t));
		offset += 4 * sizeof(int32_t);
		memcpy(target + offset, (void*)&hits_align_info.ptr->readlen, sizeof(uint32_t));
		offset += sizeof(uint32_t);
		memcpy(target + offset, (void*)&hits_align_info.ptr->score1, 4 * sizeof(uint16_t));
		offset += 4 * sizeof(uint16_t);
		memcpy(target + offset, (void*)&hits_align_info.ptr->strand, sizeof(bool));

		return target;
	} // ~Read::matchesToBinString

	// deserialize matches from string
	void unmarshallString(std::string matchStr);

	// deserialize matches from JSON
	void unmarshallJson(std::string jsonStr);
}; // ~struct Read

/**
 * Queue for Reads' records. Concurrently accessed by the Reader (producer) and the Processors (consumers)
 */
class ReadsQueue {
	int id;
	std::queue<Read> recs; // shared: Reader & Processors, Writer & Processors
	int capacity; // max size of the queue
//	int queueSizeAvr; // average size of the queue
	bool doneAdding; // flag indicating no more records will be added. Shared.
	int numPushed = 0; // shared
	int numPopped = 0; // shared
	int numPushers; // counter of threads that push reads on this queue. Used to calculate if the pushing is over.

	std::mutex lock;
	std::condition_variable cv;

public:
	ReadsQueue(int id, int capacity, int numPushers) 
		: id(id), capacity(capacity), doneAdding(false), numPushers(numPushers) {}
	~ReadsQueue() { 
		printf("Destructor called on queue %d  recs.size= %d pushed: %d  popped: %d\n", 
			id, recs.size(), numPushed, numPopped); 
	}

	void push(Read & readsrec) {
		std::unique_lock<std::mutex> l(lock);
		cv.wait(l, [this] {return recs.size() < capacity;});
		recs.push(readsrec);
		++numPushed;
//		l.unlock();
		cv.notify_one();
	}

	Read pop() {
		std::unique_lock<std::mutex> l(lock);
		cv.wait(l, [this] { return doneAdding || !recs.empty();}); //  if False - keep waiting, else - proceed.
		Read rec;
		if (!recs.empty()) {
			//		printf("%d Recs.size: %d\n", id, recs.size());
			rec = recs.front();
			recs.pop();
			++numPopped;
			if (numPopped % 10000 == 0)
				printf("\rThread %d Pushed: %d Popped: %d", std::this_thread::get_id(), numPushed, numPopped);
		}
		//	l.unlock(); // probably redundant. The lock is released when function returns
		cv.notify_one();
		return rec;
	}

	void mDoneAdding() {
		std::lock_guard<std::mutex> l(lock);
		--numPushers;
		if (numPushers == 0)
			doneAdding = true;
		cv.notify_one(); // otherwise pop can stuck not knowing the adding stopped
	}

	bool isDone() { return doneAdding && recs.empty(); }
}; // ~class ReadsQueue


// Queue for Reads IDs used to synchronize ordering of the records in Reads file and ReadStats file
struct ReadWriterCounterQueue {
	ReadWriterCounterQueue(){}
	~ReadWriterCounterQueue() {}
};

class KeyValueDatabase {
public:
	KeyValueDatabase(std::string kvdbPath) {
		// init and open key-value database for read matches
		options.IncreaseParallelism();
		options.compression = rocksdb::kXpressCompression;
		options.create_if_missing = true;
		rocksdb::Status s = rocksdb::DB::Open(options, kvdbPath, &kvdb);
		assert(s.ok());
	}
	~KeyValueDatabase() { delete kvdb; }

	void put(std::string key, std::string val)
	{
		rocksdb::Status s = kvdb->Put(rocksdb::WriteOptions(), key, val);
	}

	std::string get(std::string key)
	{
		std::string val;
		rocksdb::Status s = kvdb->Get(rocksdb::ReadOptions(), key, &val);
		return val;
	}
private:
	rocksdb::DB* kvdb;
	rocksdb::Options options;
};

// reads Reads and ReadStats files, generates Read objects and pushes them onto ReadsQueue
class Reader {
public:
	Reader(int id, ReadsQueue & readQueue, std::string & readsfile, KeyValueDatabase & kvdb, int loopCount)
		: id(id), readQueue(readQueue), readsfile(readsfile), kvdb(kvdb), loopCount(loopCount) {}
	void operator()() { read(); }
	void read();
private:
	int id;
	int loopCount; // counter of processing iterations.
	ReadsQueue & readQueue; // shared with Processor
	std::string & readsfile;
	KeyValueDatabase & kvdb; // key-value database path (from Options)
};

class Writer {
public:
	Writer(int id, ReadsQueue & writeQueue, KeyValueDatabase & kvdb)
		: id(id), writeQueue(writeQueue), kvdb(kvdb) {}
	~Writer() {}

	void operator()() { write(); }
	void write();
private:
	int id;
	ReadsQueue & writeQueue; // shared with Processor
	KeyValueDatabase & kvdb; // key-value database path (from Options)
};

class Processor {
public:
	Processor(int id, 
		ReadsQueue & readQueue,
		ReadsQueue & writeQueue,
		ReadStats & readstats,
		Index & index, 
		std::function<void(Read)> callback
	) :
		id(id), 
		readQueue(readQueue),
		writeQueue(writeQueue),
		readstats(readstats),
		index(index),
		callback(callback) {}

	void operator()() { process(); }
	void process();
private:
	int id;
	ReadsQueue & readQueue;
	ReadsQueue & writeQueue;
	ReadStats & readstats;
	Index & index;
	std::function<void(Read)> callback;
};

class Output {
public:
	Output(Runopts opts): opts(opts) {}
	~Output(){}

	void init();

private:
	Runopts & opts;
	// output streams for aligned reads (FASTA/FASTQ, SAM and BLAST-like)
	ofstream acceptedreads;
	ofstream acceptedsam;
	ofstream acceptedblast;

	char *acceptedstrings;
	char *acceptedstrings_sam;
	char *acceptedstrings_blast;
	char *logoutfile;
	char *denovo_otus_file;
	char *acceptedotumap_file;
};

// ~PARALLELTRAVERSAL_H
