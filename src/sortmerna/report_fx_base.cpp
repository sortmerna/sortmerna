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

#include <fstream>

#include "report_fx_base.h"
#include "common.hpp"
#include "options.hpp"
#include "read.hpp"
#include "readfeed.hpp"
#include "izlib.hpp"

ReportFxBase::ReportFxBase(): num_out(0), out_type(0), num_reads(0), num_hits(0), num_miss(0), num_io_bad(0), num_io_fail(0) {}

ReportFxBase::ReportFxBase(Runopts& opts): ReportFxBase()
{
	init(opts);
}

void ReportFxBase::init(Runopts& opts)
{
	validate_out_type(opts);
	set_num_out(opts);
}

void ReportFxBase::init(Readfeed& readfeed, 
                        Runopts& opts, 
                        std::vector<std::string>& fv, 
                        std::vector<std::fstream>& fsv, 
                        const std::string& fpfx, 
                        const std::string& pid_str)
{
	unsigned num_split = readfeed.num_splits * num_out;
	fsv.resize(num_split);
	fv.resize(num_split);
	// fasta/q output  WORKDIR/out/aligned_paired_fwd_0_PID.fq
	//                              pfx + sfx1 + sfx2 + sfx3 + sfx4 + ext
	INFO("num_out: ", num_out);
	for (unsigned i = 0; i < readfeed.num_splits; ++i) {
		for (unsigned j = 0, idx = 0, sense_i = 0, orig_i = 0; j < num_out; ++j) {
			std::string sfx1 = "";
			std::string sfx2 = "";

			if (num_out == 4) { // apf, apr, asf, asr
				if (j == 0) sfx1 = "_paired_fwd";
				else if (j == 1) sfx1 = "_paired_rev";
				else if (j == 2) sfx1 = "_singleton_fwd";
				else if (j == 3) sfx1 = "_singleton_rev";
			}
			else if (num_out == 2) { // ap, as | af, ar
				if (opts.is_out2) { // af, ar
					sfx2 = j == 0 ? "_fwd" : "_rev";
				}
				else if (opts.is_sout) { // ap, as
					sfx1 = j == 0 ? "_paired" : "_singleton";
				}
				else { // should never happen
					ERR("num_out = 2 implies either 'out2' or 'sout'");
					exit(1);
				}
			}

			std::string sfx3 = "_" + std::to_string(i);
			std::string sfx4 = opts.is_pid ? "_" + pid_str : "";
			std::string orig_ext = readfeed.orig_files[orig_i].isFastq ? ".fq" : ".fa";
			bool is_zip = (opts.zip_out == 1) || (readfeed.orig_files[0].isZip && opts.zip_out == -1);
			std::string gz = is_zip ? ".gz" : "";

			idx = i * num_out + j;
			fv[idx] = fpfx + sfx1 + sfx2 + sfx3 + orig_ext + gz; // e.g. aligned_paired_fwd_0_PID.fq
			if (opts.is_paired) sense_i ^= 1; // flip sense
			if (readfeed.is_two_files) orig_i ^= 1; // flip orig file index
		}
	}
}

void ReportFxBase::validate_out_type(Runopts& opts)
{
	std::stringstream ss;
	// field len      8        9        8        11        12         6    6     7      5
	std::string tl("1-file  2-files  paired  paired_in  paired_out  out2  sout  other  otype");
	std::string setf[] = { "   +    ","   +     ","   +    ","   +       ","   +        ","   +  ","   +  ","   +   " };
	std::string naf[] = { "        ","         ","        ","           ","            ","      ","      ","       " };
	std::string rules[] = {
		"no output options provided",
		"no input files provided",
		"both 1-file and 2-files specified",
		"'paired_in and paired_out cannot be used together'",
		"'sout' cannot be used with 'paired_in' or 'paired_out'"
	};

	bool is_na = false;
	int rule = -1;
	int m_1 = opts.readfiles.size() == 1 ? mask_1_file : 0x00;
	int m_2 = opts.readfiles.size() == 2 ? mask_2_file : 0x00;
	int m_p = opts.is_paired ? mask_paired : 0x00;
	int m_pin = opts.is_paired_in ? mask_paired_in : 0x00;
	int m_pout = opts.is_paired_out ? mask_paired_out : 0x00;
	int m_out2 = opts.is_out2 ? mask_out2 : 0x00;
	int m_sout = opts.is_sout ? mask_sout : 0x00;
	out_type = m_1 | m_2 | m_p | m_pin | m_pout | m_out2 | m_sout;
	if (out_type == 0x00) {
		is_na = true; rule = 0; // no options set - can never happend
	}
	if (m_1 == 0x00 && m_2 == 0x00) {
		is_na = true; rule = 1; // neither 1f or 2f specified - options validation will give an error
	}
	else if ((out_type & (mask_1_file | mask_2_file)) == (mask_1_file | mask_2_file)) {
		is_na = true; rule = 2; // both 1f and 2f specified
	}
	else if ((out_type & (mask_paired_in | mask_paired_out)) == (mask_paired_in | mask_paired_out)) {
		is_na = true; rule = 3; // both paired_in and paired_out cannot be specified
	}
	else if (out_type & mask_sout && ((out_type & mask_paired_in) || (out_type & mask_paired_out))) {
		is_na = true; rule = 4; // sout cannot be used with paired_in or paired_out
	}

	INFO("Output type: \n", tl);
	if (m_1 == mask_1_file) ss << setf[0]; else ss << naf[0];
	if (m_2 == mask_2_file) ss << setf[1]; else ss << naf[1];
	if (m_p == mask_paired) ss << setf[2]; else ss << naf[2];
	if (m_pin == mask_paired_in) ss << setf[3]; else ss << naf[3];
	if (m_pout == mask_paired_out) ss << setf[4]; else ss << naf[4];
	if (m_out2 == mask_out2) ss << setf[5]; else ss << naf[5];
	if (m_sout == mask_sout) ss << setf[6]; else ss << naf[6];
	if (opts.is_other) ss << setf[7]; else ss << naf[7];
	ss << "  " << std::hex << out_type << "\n";
	std::cout << ss.str();

	if (is_na) {
		ERR("invalid combination of output options: rule '", rule, "': '", rules[2], "' violated");
		exit(1);
	}
}

void ReportFxBase::set_num_out(Runopts& opts)
{
	if (opts.is_out2 && opts.is_sout) num_out = 4; // apf, apr, asf, asr
	else if (opts.is_out2 || opts.is_sout) num_out = 2; // ap, as | af, ar
	else num_out = 1; // a
	INFO("num_out: ", num_out);
}

/*
* write a fasta/q read
*/
void ReportFxBase::write_a_read(std::ostream& strm, Read& read, const int& dbg)
{
	std::stringstream ss;
	ss << read.header << std::endl << read.sequence << std::endl;
	if (read.format == BIO_FORMAT::FASTQ)
		ss << '+' << std::endl << read.quality << std::endl;
	if (dbg > 1) {
		try {
			strm << ss.str();
			if (strm.bad())
				num_io_bad.fetch_add(1, std::memory_order_relaxed);
			else if (strm.fail())
				num_io_fail.fetch_add(1, std::memory_order_relaxed);
			else if (strm.good())
				num_reads.fetch_add(1, std::memory_order_relaxed);

			if (read.is_hit)
				num_hits.fetch_add(1, std::memory_order_relaxed);
			else
				num_miss.fetch_add(1, std::memory_order_relaxed);
		}
		catch (const std::exception& e) {
			ERR("failed writing to stream. Num reads processed so far: ", num_reads, " Current read id: ", read.id, " - ", e.what());
			exit(1);
		}
	}
	else {
		strm << ss.str();
	}
}


void ReportFxBase::write_a_read(std::ostream& strm, Read& read, Readstate& rstate, Izlib& izlib, bool is_last, const int& dbg)
{
	++rstate.read_count;
	std::stringstream ss;
	if (is_last && read.sequence.empty()) {
		ss.str("");
	}
	else {
		ss << read.header << std::endl << read.sequence << std::endl;
		if (read.format == BIO_FORMAT::FASTQ)
			ss << '+' << std::endl << read.quality << std::endl;
	}
	auto ret = izlib.defstr(ss.str(), strm, is_last); // Z_STREAM_END | Z_OK - ok
	if (ret < Z_OK || ret > Z_STREAM_END) {
		ERR("Failed deflating readstring: ", ss.str(), " zlib status: ", ret);
	}
}