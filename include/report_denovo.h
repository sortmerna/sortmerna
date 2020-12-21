#pragma once

#include <string>
#include "report.h"

// forward
class Read;

class ReportDenovo : public Report
{
	std::string ext = ".fa";
public:
	ReportDenovo(Runopts& opts);
	ReportDenovo(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void append(int id, Read& read, Runopts& opts);
};