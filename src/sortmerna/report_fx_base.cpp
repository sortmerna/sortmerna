#include <fstream>

#include "report_fx_base.h"
#include "common.hpp"
#include "options.hpp"
#include "read.hpp"
#include "readfeed.hpp"
#include "izlib.hpp"

ReportFxBase::ReportFxBase(): num_out(0), out_type(0) {}

ReportFxBase::ReportFxBase(Runopts& opts): ReportFxBase()
{
	init(opts);
}

void ReportFxBase::init(Runopts& opts)
{
	calc_out_type(opts);
	set_num_out();
}

void ReportFxBase::init(Readfeed& readfeed, Runopts& opts, std::vector<std::string>& fv, std::vector<std::fstream>& fsv, std::string& fpfx, std::string& pid_str)
{
	auto num_split = readfeed.num_splits * num_out;
	fsv.resize(num_split);
	fv.resize(num_split);
	// fasta/q output  WORKDIR/out/aligned_paired_fwd_0_PID.fq
	//                              pfx + sfx1 + sfx2 + sfx3 + sfx4 + ext
	for (int i = 0; i < readfeed.num_splits; ++i) {
		for (int j = 0, idx = 0, orig_idx = 0; j < num_out; ++j) {
			std::string sfx1 = "";
			if (out_type == 0x44 || out_type == 0x46) { // apf, apr, asf, asr
				if (j == 0 || j == 1) sfx1 = "_paired";
				else if (j == 2 || j == 3) sfx1 = "_singleton";
			}
			//(out_type == 0x5C || out_type == 0x5E) { // af, ar, of, or
			std::string sfx2 = "";
			sfx2 = opts.is_out2 && orig_idx == 0 ? "_fwd" : "_rev";

			std::string sfx3 = "_" + std::to_string(i);
			std::string sfx4 = opts.is_pid ? "_" + pid_str : "";
			std::string orig_ext = readfeed.orig_files[orig_idx].isFastq ? ".fq" : ".fa";

			// test the file(s)
			idx = i * num_out + j;
			fv[idx] = fpfx + sfx1 + sfx2 + sfx3 + orig_ext; // e.g. aligned_paired_fwd_0_PID.fq
			orig_idx = readfeed.is_two_files ? orig_idx ^= 1 : orig_idx; // flip
		}
	}
}

void ReportFxBase::calc_out_type(Runopts& opts)
{
	if (opts.readfiles.size() == 1) out_type |= mask_1_file;
	if (opts.is_paired) out_type |= mask_paired;
	if (opts.readfiles.size() == 2) out_type |= mask_2_file;
	if (opts.is_other) out_type |= mask_other;
	if (opts.is_paired_in) out_type |= mask_paired_in;
	if (opts.is_paired_out) out_type |= mask_paired_out;
	if (opts.is_out2) out_type |= mask_out2;

	INFO("Output type: ", out_type);
}

void ReportFxBase::set_num_out()
{
	if (out_type == 0x04 || out_type == 0x06) num_out = 2; // ap, as
	else if (out_type == 0x14 || out_type == 0x16
		|| out_type == 0x24 || out_type == 0x26) num_out = 1; // a
	else if (out_type == 0x44 || out_type == 0x46) num_out = 4; // apf, apr, asf, asr   'af.size != ar.size' in general. Only aligned reads.
	else if (out_type == 0x4C || out_type == 0x4E) num_out = 4; // 8 if other
	else if (out_type == 0x5C || out_type == 0x5E) num_out = 2; // 4 if other: af, ar, of, or
}

/*
* write a fasta/q read
*/
void ReportFxBase::write_a_read(std::ostream& strm, Read& read)
{
	strm << read.header << std::endl << read.sequence << std::endl;
	if (read.format == BIO_FORMAT::FASTQ)
		strm << '+' << std::endl << read.quality << std::endl;
}


void ReportFxBase::write_a_read(std::ostream& strm, Read& read, Readstate& rstate, Izlib& izlib)
{
	++rstate.read_count;
	std::stringstream ss;
	ss << read.header << std::endl << read.sequence << std::endl;
	auto ret = izlib.defstr(ss.str(), strm); // Z_STREAM_END | Z_OK - ok
	if (ret < Z_OK || ret > Z_STREAM_END) {
		ERR("Failed deflating readstring: ", ss.str(), " zlib status: ", ret);
	}
}