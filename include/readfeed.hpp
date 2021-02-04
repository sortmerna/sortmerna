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
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
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
	void init(std::vector<std::string>& readfiles);
	void init_split_files();
	void init_reading();
	void init_vzlib_in();
	//void init_vstate_in();
	bool split();
	bool define_format();
	void count_reads();
	void write_descriptor();
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
	unsigned num_sense; // number of read's senses (fwd/rev)
	unsigned num_reads_tot; // count of reads in all streams
	unsigned length_all; // length of all reads from all files
	unsigned min_read_len;
	unsigned max_read_len;
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