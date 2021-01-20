#include "report_denovo.h"
#include "common.hpp"
#include "options.hpp"
#include "read.hpp"
#include "references.hpp"
#include "refstats.hpp"
#include "readfeed.hpp"

ReportDenovo::ReportDenovo(Runopts& opts) : Report(opts), base() {}

ReportDenovo::ReportDenovo(Readfeed& readfeed, Runopts& opts) : ReportDenovo(opts)
{
	init(readfeed, opts);
}

void ReportDenovo::init(Readfeed& readfeed, Runopts& opts)
{
	base.init(opts);
	base.init(readfeed, opts, fv, fsv, opts.aligned_pfx.string() + "_denovo", pid_str);
	openfw(); // open output files for writing
	is_zip = readfeed.orig_files[0].isZip;
	if (is_zip)	init_zip();

	//fv.resize(readfeed.num_splits);
	//fsv.resize(readfeed.num_splits);
	//is_zip = readfeed.orig_files[0].isZip;
	// WORKDIR/out/aligned_denovo_0_PID.fa
	//for (int i = 0; i < readfeed.num_splits; ++i) {
	//	std::string sfx1 = "_" + std::to_string(i);
	//	std::string sfx2 = opts.is_pid ? "_" + pid_str : "";
	//	std::string gz = is_zip ? ".gz" : "";
	//	fv[i] = opts.aligned_pfx.string() + "_denovo" + sfx1 + sfx2 + ext + gz;
	//	openfw(i);
	//}
	//if (is_zip) init_zip();
}

void ReportDenovo::append(int id, std::vector<Read>& reads, Runopts& opts, bool is_last)
{
	if (opts.is_paired) {
		if (!reads[0].is_denovo && !reads[1].is_denovo)
			return; // neiher is denovo - nothing to write

		// caclulate the index of the output file to write to
		for (int i = 0, idx = 0; i < reads.size(); ++i)
		{
			// 1 output file a (aligned reads)
			if (base.num_out == 1) {
				if (reads[i].is_denovo || opts.is_paired_in)
					idx = id;
				else
					continue;
			}
			// 2 output files dp,ds (sout) | df,dr (out2)
			else if (base.num_out == 2) {
				if (opts.is_out2) {
					if (opts.is_paired_out && !(reads[0].is_denovo && reads[1].is_denovo))
						break; // if not both aligned -> non-aligned
					else if (opts.is_paired_in || reads[i].is_denovo)
						idx = id * base.num_out + i;
				}
				else if (opts.is_sout) {
					if (reads[0].is_denovo && reads[1].is_denovo)
						idx = id * base.num_out; // both to 'dp' [0]
					else if (reads[i].is_denovo)
						idx = id * base.num_out + 1; // hit to 'ds' [1]
					else
						continue; // ignore non-aligned read
				}
			}
			// 4 output files: dpf, dpr, dsf, dsr
			// being here means both is_out2 & is_sout were set => No paired_in/out
			else if (base.num_out == 4) {
				if (reads[0].is_denovo && reads[1].is_denovo)
					idx = id * base.num_out + i; // 0 -> apf, 1 -> apr
				else if (reads[i].is_denovo)
					idx = id * base.num_out + i + 2; // 0 -> asf, 1 -> asr
				else
					continue; // ignore a non-aligned singleton
			}
			else {
				ERR("min number of output files is 1, max number of output files is 4. The current value is ", base.num_out);
				exit(1);
			}

			if (is_zip)
				base.write_a_read(fsv[idx], reads[i], vstate_out[idx], vzlib_out[idx], is_last);
			else
				base.write_a_read(fsv[idx], reads[i]);
		}
	}//~if paired
	// non-paired
	else
	{
		// the read was accepted - output
		if (reads[0].is_denovo)
		{
			if (is_zip)
				base.write_a_read(fsv[0], reads[0], vstate_out[0], vzlib_out[0], is_last);
			else
				base.write_a_read(fsv[0], reads[0]);
		}
	}
}

void ReportDenovo::merge(int num_splits)
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
		strip_path_sfx(fv[i]);
	}
}