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
			   Laurent No�      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mika�l Salson    mikael.salson@lifl.fr
			   H�l�ne Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

/*
 * FILE: readfeed.hpp
 * Created: Nov 06, 2017 Mon
 *
 * Encapsulates operations on the 'Reads' files both flat and archived
 */

#pragma once

#include <string>
#include <fstream> // std::ifstream
#include <filesystem>

#include "common.hpp"
#include "izlib.hpp"
#include "readstate.h"
#include "readfile.h"

 // forward
class Read;
class KeyValueDatabase;
struct Runopts;

/* 
 * Database-like interface for accessing reads
 */
class Readfeed {
public:
	Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, std::filesystem::path& basedir, bool is_paired);
	Readfeed(FEED_TYPE type, std::vector<std::string>& readfiles, const unsigned num_parts, std::filesystem::path& basedir, bool is_paired);

	void run();
	bool next(int inext, std::string& readstr);
	void reset();
	void rewind();
	void rewind_in();
	void init(std::vector<std::string>& readfiles, const int& dbg = 0);
	void init_split_files();
	void init_reading();
	void init_vzlib_in();
	//void init_vstate_in();
	bool split();
	/*
     * - define input files format (FASTA, FASTQ) and compression (gz, non-gz)
     * - test reading by getting one read
     *
     * logic:
     *   gz:
   	 *    1F 8B 08 08 61 78 C5 5E 00 03 53 52 52 31 36 33 35 38 36 34 5F 31 5F 35 4B 2E 66 61 73 74 71 00
   	 *     0 byte at 8th position_|  |  |_file name starts in ASCII                   file name end_|  |_ 0 byte
   	 *                      ETX byte_|
   	 * flat:
   	 *   first 100 bytes are ascii (<=127 x7F), first char is '@' (x40), and can infer fasta or fastq
    */
	bool define_format(const int& dbg = 0);
	void count_reads();
	void write_descriptor();
	/*
     * verify the split was already performed and the feed is ready
     * logic:
	 *   Read feed descriptor and verify:
	 *     - original files have the same name, size and line count
	 *     - num splits are the same
	 *     - split files names, count, and line count are the same
     * Split readfeed descriptor:
	 *   timestamp: xxx
	 *   num_input: 2  # number of input files
	 *   num_parts: 3  # number of split parts. split[].size = num_input * num_parts
	 *   input:
	 *     - file_1: name, sha
	 *     - file_2: name, sha
	 * split:
	 *  - file_1: name, sha
	 *  - file_2: name, sha
	 *  ...
	 *  - file_n: name, sha
    */
	bool is_split_ready();
	/*
	* delete the readfeed i.e. delete Split files (Not original reads files) described in the readfeed descriptor
	* @return  count of deleted split files
	*/
	int clean();
	static bool hasnext(std::ifstream& ifs);
	static bool loadReadByIdx(Read& read);
	static bool loadReadById(Read& read);

private:
	/*
     * get a next read string from a reads file
     * 20201012: TODO: exactly the same as next(int, uint, str, bool) but doesn't count sequence length.
	 *  			   Counting sequence length could be a very small overhead -> no need for this overload?
     * 
     * @param inext   IN  index of the stream to read.
     * @param readstr IN  read sequence
     * @param is_orig IN  flags to return the original read string. If false, then the read string has format: 'read_id \n header \n sequence [\n quality]'
     * @param files   IN  array of read files' descriptors
     * @return            true if record exists, else false
    */
	bool next(int inext, std::string& readstr, bool is_orig, std::vector<Readfile>& files);
	//bool next(int inext, std::string& readstr, unsigned& readlen, bool is_orig = false); // \n separated read data. FA - 2 lines, FQ - 4 lines
	/*
     * uses zlib if files are gzipped or reads flat files
     *
     * @param inext   IN  index into the Readfeed::files vector
	 * @param readstr IN  read string
     * @param readlen IN  length of the read sequence
     * @param is_orig IN  flags to return the original read string. If false, then the read string has format: 'read_id \n header \n sequence [\n quality]'
     * @param files   IN  array of read files' descriptors
	 * @return            true if record exists, else false
    */
	bool next(int inext, std::string& readstr, unsigned& readlen, bool is_orig, std::vector<Readfile>& files);

public:
	FEED_TYPE type;
	bool is_done; // flags end of all read streams
	bool is_ready; // flags the read feed is ready i.e. no need to run split
	bool is_format_defined; // flags the file format is defined i.e. 'define_format' was success
	bool is_two_files; // flags two read files are processed (otherwise single file)
	bool is_paired;
	unsigned num_orig_files; // number of original reads files
	unsigned num_splits;
	unsigned num_split_files;
	uint32_t num_sense; // number of read's senses (fwd/rev) i.e. max 2
	uint64_t num_reads_tot; // count of reads in all streams
	uint64_t length_all; // length of all reads from all files
	uint32_t min_read_len;
	uint32_t max_read_len;
	std::filesystem::path& basedir; // root directory for split files (opts.readb)
	std::vector<Readfile> orig_files;
private:
	std::vector<Readfile> split_files;

	// input processing
	std::vector<std::ifstream> ifsv; // [fwd_0, rev_0, fwd_1, rev_1, ... fwd_n-1, rev_n-1]
	std::vector<Izlib> vzlib_in;
	std::vector<Readstate> vstate_in;

	// output processing
	std::vector<std::ofstream> ofsv;
	std::vector<Izlib> vzlib_out;
	std::vector<Readstate> vstate_out;
};

// ~readfeed.hpp