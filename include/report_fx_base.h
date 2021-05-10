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

#pragma once

#include <fstream>
#include <vector>
#include <atomic>

// forward
class Read;
struct Readstate;
class Readfeed;
struct Runopts;
class Izlib;

/*
* FASTX report's data and functionality common to both FASTX aligned and other (non-aligned) reports
*/
class ReportFxBase {
public:
	ReportFxBase();
	ReportFxBase(Runopts& opts);
	void init(Runopts& opts);
	/*
	* init output file names. It has a different semantics than init(opts) above.
	*/
	void init(Readfeed& readfeed, Runopts& opts, std::vector<std::string>& fv, std::vector<std::fstream>& fsv, const std::string& fpfx, const std::string& pid_str);
	void write_a_read(std::ostream& strm, Read& read, const int& dbg=0);
	void write_a_read(std::ostream& strm, Read& read, Readstate& rstate, Izlib& izlib, bool is_last=false, const int& dbg = 0);

	unsigned num_out; // number of aligned output files (1 | 2 | 4) depending on the output type below
	/*
	* output type: 255/2 -> ~120 possible types
	* defines the number of the files to output
	* out_type = mask | mask | mask ... where 'mask' depends on the specified option
	*/
	int out_type;
	std::atomic<unsigned long> num_reads; // count of reads processed so far. Use For debugging.
	std::atomic<unsigned long> num_hits; // count of hits
	std::atomic<unsigned long> num_miss; // count of misses (other)
	std::atomic<unsigned long> num_io_bad;
	std::atomic<unsigned long> num_io_fail;

private:
	/*
	* 1-file  2-files  paired  paired_in  paired_out  out2  sout  other
	* x01     x02      x04     x08        x10         x20   x40   x80
	*/
	const int mask_1_file = 0x01;
	const int mask_2_file = 0x02;
	const int mask_paired = 0x04;
	const int mask_paired_in = 0x08;
	const int mask_paired_out = 0x10;
	const int mask_out2 = 0x20;
	const int mask_sout = 0x40;
	const int mask_other = 0x80;

	/*
	* validate the output type
	*/
	void validate_out_type(Runopts& opts);
	/*
	* set number of output files: 1 | 2 | 4 - calculated depending on the selected options
	*/
	void set_num_out(Runopts& opts);
};
