#include "report_denovo.h"
#include "common.hpp"
#include "options.hpp"
#include "read.hpp"
#include "references.hpp"
#include "refstats.hpp"
#include "readfeed.hpp"

ReportDenovo::ReportDenovo(Runopts& opts) : Report(opts) {}

ReportDenovo::ReportDenovo(Readfeed& readfeed, Runopts& opts) : ReportDenovo(opts)
{
	init(readfeed, opts);
}

void ReportDenovo::init(Readfeed& readfeed, Runopts& opts)
{
	fv.resize(readfeed.num_splits);
	fsv.resize(readfeed.num_splits);
	is_zip = readfeed.orig_files[0].isZip;
	// WORKDIR/out/aligned_denovo_0_PID.fa
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + pid_str : "";
		std::string gz = is_zip ? ".gz" : "";
		fv[i] = opts.aligned_pfx.string() + "_denovo" + sfx1 + sfx2 + ext + gz;
		openfw(i);
	}
	if (is_zip) init_zip();
}

void ReportDenovo::append(int id, Read& read, Runopts& opts)
{
	std::stringstream ss;
	ss << read.header << std::endl << read.sequence << std::endl;
	if (is_zip) {
		auto ret = vzlib_out[id].defstr(ss.str(), fsv[id]); // Z_STREAM_END | Z_OK - ok
		if (ret < Z_OK || ret > Z_STREAM_END) {
			ERR("Failed deflating readstring: ", ss.str(), " zlib status: ", ret);
		}
	}
	else
		fsv[id] << ss.str();
}