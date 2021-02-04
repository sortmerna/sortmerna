/*
 @copyright 2016-2021  Clarity Genomics BVBA
 @copyright 2012-2016  Bonsai Bioinformatics Research Group
 @copyright 2014-2016  Knight Lab, Department of Pediatrics, UCSD, La Jolla

 @parblock
 SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 This is a free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SortMeRNA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SortMeRNA. If not, see <http://www.gnu.org/licenses/>.
 @endparblock

 @contributors Jenya Kopylova   jenya.kopylov@gmail.com
			   Laurent Noé      laurent.noe@lifl.fr
			   Pierre Pericard  pierre.pericard@lifl.fr
			   Daniel McDonald  wasade@gmail.com
			   Mikaël Salson    mikael.salson@lifl.fr
			   Hélène Touzet    helene.touzet@lifl.fr
			   Rob Knight       robknight@ucsd.edu
*/

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