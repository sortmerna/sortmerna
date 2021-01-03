#pragma once
#include "report.h"
#include "report_fx_base.h"

// forward
class Readfead;
class ReportFxBase;

class ReportFxOther : public Report {
public:
	ReportFxOther(Runopts& opts);
	ReportFxOther(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void merge(int num_splits) override;
	void append(int id, std::vector<Read>& reads, Runopts& opts, bool is_last=false);
	ReportFxBase& getBase();

private:
	ReportFxBase base;
};