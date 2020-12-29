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
	* @param id  processing thread enumerator 0,1,2..
    * @param reads: 1 or 2 (paired) reads
	* @param opts
	* @param is_last   flags the passed read is the last
    */
	void append(int id, std::vector<Read>& reads, Runopts& opts, bool is_last=false);
	ReportFxBase& getBase();

private:
	ReportFxBase base;
};
