#pragma once

#include <string>
#include "report.h"

// forward
class Read;
class References;

class ReportBiom : public Report
{
	std::string ext = ".sam";
public:
	ReportBiom(Runopts& opts);
	ReportBiom(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void append();
};
