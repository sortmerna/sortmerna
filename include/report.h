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
	* blast, sam, denovo - use base implementation
	* fastx - override
	*/
	virtual void merge(int num_splits);
	/*
	* init zip interface. So far common to all the reports.
	*/
	virtual void init_zip();

	void openfr(int idx); // open a single file for reading
	void openfw(int idx);
	void openfw(); // open report files for writing
	void closef(int idx); // close a file
	void closef(); // close report files
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
