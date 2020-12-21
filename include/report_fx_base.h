#pragma once

#include <fstream>
#include <vector>

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
	void init(Readfeed& readfeed, Runopts& opts, std::vector<std::string>& fv, std::vector<std::fstream>& fsv, std::string& fpfx, std::string& pid_str);
	void write_a_read(std::ostream& strm, Read& read);
	void write_a_read(std::ostream& strm, Read& read, Readstate& rstate, Izlib& izlib);

	int num_out; // number of aligned output files
	int out_type; //  = 0x00

private:
	const int mask_1_file = 0x01;
	const int mask_paired = 0x02;
	const int mask_2_file = 0x04;
	const int mask_other = 0x08;
	const int mask_paired_in = 0x10;
	const int mask_paired_out = 0x20;
	const int mask_out2 = 0x40;

	void calc_out_type(Runopts& opts);
	void set_num_out();
};
