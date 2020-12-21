#pragma once
#include "report.h"
#include "report_fx_base.h"

// forward
class Read;

class ReportFastx: public Report
{
public:
	ReportFastx(Runopts& opts);
	ReportFastx(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void merge(int num_splits) override;
	/*
    * called on each read or each pair of reads (if paired)
    * writes to aliged.fasta
    *
    * @param reads: 1 or 2 (paired) reads
    */
	void append(int id, std::vector<Read>& reads, Runopts& opts);
	ReportFxBase& getBase();

private:
	ReportFxBase base;
};
