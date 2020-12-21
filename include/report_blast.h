#pragma once
#include <string>
#include "report.h"

// forward
class References;
class Refstats;
class Read;

class ReportBlast : public Report
{
	std::string ext = ".blast";
public:
	ReportBlast(Runopts& opts);
	ReportBlast(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void append(int id, Read& read, References& refs, Refstats& refstats, Runopts& opts);
};