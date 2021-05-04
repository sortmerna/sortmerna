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

#include <vector>
#include <fstream>

#include "readstate.h"
#include "readfile.h"
#include "izlib.hpp"

// forward
class Readfeed;
struct Runopts;

/*
 * base class for writing reports
*/
class Report
{
public:
	Report(Runopts& opts);
	virtual ~Report() = 0;
	virtual void init(Readfeed& readfeed, Runopts& opts) = 0;
	/*
	* merge split report files
	* override if necessary
	*/
	//virtual void merge(int num_splits, const int& dbg=0);
	virtual void merge(int num_splits, const int& num_out, const int& dbg = 0);
	/*
	* init zip interface. So far common to all the reports.
	*/
	virtual void init_zip();

	void openfr(unsigned idx); // open a single file for reading
	/*
	* open report files for writing
	* @param dbg  debug level, see 'Runopts.dbg_level'
	*/
	void openfw(const int& dbg=0);
	void openfw(size_t idx, const int& dbg=0);
	void closef(const int& dbg=0); // close report files
	void closef(size_t idx, const int& dbg=0); // close a file
	int finish_deflate();
	/*
	* strip the split suffix from the output file name 
	* and rename the final output e.g. 
	*   'aligned_0.blast -> aligned.blast'
	* @param path  e.g. $HOME/sortmerna/run/out/aligned_0.blast
	* @param sfx   e.g. '_0'
	*/
	void strip_path_sfx(std::string& path, std::string sfx="_0");

protected:
	std::string pid_str; // std::to_string(getpid());
	bool is_zip; // flags the report is compressed
	std::vector<std::string> fv; // report files
	std::vector<std::fstream> fsv;

	// reading compressed out files when merging the final output
	std::vector<Izlib> vzlib_in;
	std::vector<Readstate> vstate_in;

	// writing compressed output during report generation and the final compressed output
	std::vector<Izlib> vzlib_out;
	std::vector<Readstate> vstate_out;
};
