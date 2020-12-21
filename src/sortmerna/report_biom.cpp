#include "report_biom.h"
#include "common.hpp"

ReportBiom::ReportBiom(Runopts& opts) : Report(opts) {}

ReportBiom::ReportBiom(Readfeed& readfeed, Runopts& opts) : ReportBiom(opts)
{
	init(readfeed, opts);
}

void ReportBiom::init(Readfeed& readfeed, Runopts& opts)
{
	INFO("TODO");
}

void ReportBiom::append()
{
	openfw(0);
	fsv[0] << "\"id:\"null,";
	fsv[0] << "\"format\": \"Biological Observation Matrix 1.0.0\",";
	fsv[0] << "\"format_url\": \"http://biom-format.org/documentation/format_versions/biom-1.0.html\"";
	fsv[0] << "\"type\": \"OTU table\",";
	fsv[0] << "\"generated_by\": \"SortMeRNA v2.0\",";
	fsv[0] << "\"date\": \"\",";
	fsv[0] << "\"rows\":[";
	fsv[0] << "\"matrix_type\": \"sparse\",";
	fsv[0] << "\"matrix_element_type\": \"int\",";
	fsv[0] << "\"shape\":";
	fsv[0] << "\"data\":";
	fsv[0].close();
}