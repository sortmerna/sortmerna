#include <filesystem>

#include "common.hpp"
#include "report_fastx.h"
#include "readfeed.hpp"
#include "options.hpp"
#include "read.hpp"

ReportFastx::ReportFastx(Runopts& opts): Report(opts), base(ReportFxBase()) {}
ReportFastx::ReportFastx(Readfeed& readfeed, Runopts& opts): ReportFastx(opts)
{ 
	init(readfeed, opts); 
}

void ReportFastx::init(Readfeed& readfeed, Runopts& opts)
{
	base.init(opts);
	base.init(readfeed, opts, fv, fsv, opts.aligned_pfx.string(), pid_str);
	openfw(); // open output files for writing

	// prepare zlib interface for writing split files
	vzlib_out.resize(fv.size(), Izlib(true, true));
	for (auto i = 0; i < vzlib_out.size(); ++i) {
		vzlib_out[i].init(true);
	}

	// prepare Readstates OUT
	vstate_out.resize(fv.size());
} // ~ReportFastx::init

void ReportFastx::append(int id, std::vector<Read>& reads, Runopts& opts)
{
	if (reads.size() == 2) {
		if (opts.is_paired_in) {
			// if Either is aligned -> both to aligned
			if (reads[0].is_hit || reads[1].is_hit) {
				// validate the reads are paired in case of two reads files
				if (opts.readfiles.size() == 2
					&& (reads[0].read_num != reads[1].read_num
						|| reads[0].readfile_idx == reads[1].readfile_idx))
				{
					ERR("Paired validation failed: reads[0].id= ", reads[0].id, " reads[0].read_num = ",
						reads[0].read_num, " reads[0].readfile_idx= ", reads[0].readfile_idx,
						" reads[1].id=", reads[1].id, " reads[1].read_num = ", reads[1].read_num,
						" reads[1].readfile_idx = ", reads[1].readfile_idx);
					exit(EXIT_FAILURE);
				}

				// reads[0]
				for (int i = 0, idx = id * reads.size(); i < reads.size(); ++i)
				{
					// fwd and rev go into the same file if not OUT2 else to different files
					base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx]);
					if (opts.is_out2) {
						++idx; // fwd and rev go into different files
					}
				}
			}
		}
		else if (opts.is_paired_out) {
			// if Both aligned -> aligned
			if (reads[0].is_hit && reads[1].is_hit) {
				for (int i = 0, idx = id * reads.size(); i < reads.size(); ++i)
				{
					if (opts.is_out2) {
						base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx]); // fwd and rev go into different files
						++idx;
					}
					else {
						base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx]); // fwd and rev go into the same file
					}
				}
			}
		}
		else {
			// Neither 'paired_in' nor 'paired_out' specified:
			//   aligned     -> aligned file
			//   non-aligned -> other file
			auto idx = id * base.num_out;
			if (reads[0].is_hit && reads[1].is_hit) {
				if (opts.is_out2) {
					auto ii = (size_t)idx + 1;
					base.write_a_read(fsv[idx], reads[0], vstate_out[idx], vzlib_out[idx]); // apf
					base.write_a_read(fsv[ii], reads[1], vstate_out[ii], vzlib_out[ii]); // apr
				}
				else {
					base.write_a_read(fsv[idx], reads[0], vstate_out[idx], vzlib_out[idx]); // ap
					base.write_a_read(fsv[idx], reads[1], vstate_out[idx], vzlib_out[idx]); // ap
				}
			}
			else if (reads[0].is_hit) { // fwd hit
				if (opts.is_out2) {
					auto ii = (size_t)idx + 2;
					base.write_a_read(fsv[ii], reads[0], vstate_out[ii], vzlib_out[ii]); // asf
				}
				else {
					base.write_a_read(fsv[idx], reads[0], vstate_out[idx], vzlib_out[idx]); // as
				}
			}
			else if (reads[1].is_hit) { // rev hit
				if (opts.is_out2) {
					auto ii = (size_t)idx + 3;
					base.write_a_read(fsv[ii], reads[0], vstate_out[ii], vzlib_out[ii]); // asr
				}
				else {
					auto ii = (size_t)idx + 1;
					base.write_a_read(fsv[ii], reads[0], vstate_out[ii], vzlib_out[ii]); // as
				}
			}
		}
	}//~if paired
	// non-paired
	else
	{
		// the read was accepted - output
		if (reads[0].is_hit)
		{
			base.write_a_read(fsv[0], reads[0], vstate_out[0], vzlib_out[0]);
		}
	}
} // ~ReportFasta::append

void ReportFastx::merge(int num_splits)
{
	for (int i = 0; i < base.num_out; ++i) {
		openfw(i);
		for (int j = 1; j < num_splits; ++j) {
			auto idx = i + j * base.num_out;
			openfr(idx);
			fsv[i] << fsv[idx].rdbuf();
			INFO("merged ", fv[idx], " -> ", fv[i]);
			closef(idx);
			std::filesystem::remove(fv[idx]);
			INFO("deleted ", fv[idx]);
		}
		closef(i);
	}
}

ReportFxBase& ReportFastx::getBase()
{
	return base;
}