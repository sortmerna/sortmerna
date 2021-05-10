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

#pragma once

#include <string>
#include "report.h"
#include <atomic>

// forward
class References;
class Refstats;
class Read;

class ReportBlast : public Report
{
	std::string ext = ".blast";

public:
	// debug
	std::atomic<uint64_t> n_aligned; // reads passing E-value threshold.
	std::atomic<uint64_t> n_yid_ncov; // SW + ID - COV i.e. aligned passing ID, failing COV
	std::atomic<uint64_t> n_nid_ycov; // SW - ID + COV i.e. aligned failing ID, passing COV
	std::atomic<uint64_t> n_yid_ycov; // [2] SW + ID + COV i.e. aligned passing ID, passing COV
	std::atomic<uint64_t> n_denovo; // [4] SW - ID - COV i.e. 'de novo' reads, aligned failing ID, failing COV

	ReportBlast(Runopts& opts);
	ReportBlast(Readfeed& readfeed, Runopts& opts);
	void init(Readfeed& readfeed, Runopts& opts) override;
	void append(const uint32_t& id, Read& read, References& refs, Refstats& refstats, Runopts& opts);
};