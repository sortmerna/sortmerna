#pragma once

#include <string>
#include "report.h"

// forward
class Read;
class References;

class ReportSam : public Report
{
	std::string ext = ".sam";
public:
	ReportSam(Runopts& opts);
	ReportSam(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void append(int id, Read& read, References& refs, Runopts& opts);
	void write_header(Runopts& opts);
};
