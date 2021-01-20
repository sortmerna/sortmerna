#pragma once

#include "report.h"
#include "report_fx_base.h"

// forward
class Read;

class ReportDenovo : public Report
{
	std::string ext = ".fa";
public:
	ReportDenovo(Runopts& opts);
	ReportDenovo(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void merge(int num_splits) override;
	void append(int id, std::vector<Read>& reads, Runopts& opts, bool is_last = false);

private:
	ReportFxBase base;
};