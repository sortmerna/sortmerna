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
	// WORKDIR/out/aligned_denovo_0_PID.fa
	for (int i = 0; i < readfeed.num_splits; ++i) {
		std::string sfx1 = "_" + std::to_string(i);
		std::string sfx2 = opts.is_pid ? "_" + pid_str : "";
		fv[i] = opts.aligned_pfx.string() + "_denovo" + sfx1 + sfx2 + ext;
		INFO("Testing file: ", std::filesystem::absolute(std::filesystem::path(fv[i])));
		fsv[i].open(fv[i]);
		fsv[i].close();
		if (!fsv[i]) {
			ERR("Failed stream on file ", fv[i]);
			exit(EXIT_FAILURE);
		}
	}
}

void ReportDenovo::append(int id, Read& read, Runopts& opts)
{
	fsv[id] << read.header << std::endl << read.sequence << std::endl;
}