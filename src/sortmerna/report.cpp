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

#include <string>
#include <filesystem>

#include "unistd.h"
#include "report.h"
#include "common.hpp"
#include "options.hpp"

Report::Report(Runopts& opts) : pid_str(std::to_string(getpid())), is_zip(false) {}
Report::~Report() {	closef(); }

void Report::init_zip()
{
	// prepare zlib interface for writing split files
	vzlib_out.resize(fv.size(), Izlib(true, true));
	for (auto& zlibm: vzlib_out) {
		zlibm.init(true);
	}

	// prepare Readstates OUT
	vstate_out.resize(fv.size());
}

void Report::merge(const uint32_t& num_splits, const uint32_t& num_out, const int& dbg)
{
	for (uint32_t i = 0; i < num_out; ++i) {
		if (!fsv[i].is_open()) { 
			openfw(i, dbg);
		}
		else {
			fsv[i].seekp(0); //rewind
		}

		for (uint32_t j = 1; j < num_splits; ++j) {
			uint32_t idx = i + j * num_out;
			auto fsz = std::filesystem::file_size(fv[idx]);
			if (dbg > 1) {
				INFO("input file idx: ", idx, " size: ", fsz);
				if (fsz == 0)
					INFO("skipping empty file at idx: ", idx);
			}
			if (fsz > 0) {
				if (!fsv[idx].is_open()) {
					openfr(idx);
				}
				else {
					fsv[idx].seekg(0);
				}
				if (dbg > 1) {
					INFO("merging ifs@idx: ", idx);
				}
				fsv[i] << fsv[idx].rdbuf();
				INFO("merged ", fv[idx], " -> ", fv[i]);
				if (dbg > 1) {
					INFO("merged ifs@idx: ", idx, " ifs bad: ", fsv[idx].bad(), " end ifs pos: ", fsv[idx].tellg(), " to ofs@idx: ", i, " ofs bad: ", fsv[i].bad());
				}
				fsv[idx].close();
			}
			std::filesystem::remove(fv[idx]);
			INFO("deleted ", fv[idx]);
		}
		fsv[i].close();
		strip_path_sfx(fv[i]);
	}
} // ~Report::merge

void Report::openfw(const uint32_t& idx, const int& dbg)
{
	if (!fsv[idx].is_open()) {
		fsv[idx].open(fv[idx], std::ios::binary | std::ios::app);
	}
	if (!fsv[idx].good()) {
		ERR("Could not open output file number [", idx, "] : [", fv[idx], "] for writing.");
		exit(EXIT_FAILURE);
	}
	else {
		if (dbg > 0)
			INFO("Opened output file ", fv[idx], " for writing.");
	}
}

void Report::openfw(const int& dbg)
{
	for (size_t i = 0; i < fv.size(); ++i) {
		openfw(i, dbg);
	}
}

void Report::openfr(unsigned idx)
{
	if (!fsv[idx].is_open()) {
		fsv[idx].open(fv[idx], std::ios::binary | std::ios::in);
	}
	if (!fsv[idx].good()) {
		ERR("Could not open output file [", fv[idx], "] for reading.");
		exit(EXIT_FAILURE);
	}
	else {
		INFO("Opened output file ", fv[idx], " for reading.");
	}
}

void Report::closef(const uint32_t& idx, const int& dbg)
{
	if (fsv[idx].is_open()) {
		fsv[idx].flush();
		fsv[idx].close();
		if (dbg > 1)
			INFO("Closed output file: ", fv[idx]);
	}
}

void Report::closef(const int& dbg)
{
	for (uint32_t i = 0; i < fsv.size(); ++i) {
		closef(i, dbg);
	}
}

int Report::finish_deflate()
{
	int ret = 0;
	if (is_zip) {
		for (unsigned i = 0; i < vzlib_out.size(); ++i) {
			ret += vzlib_out[i].finish_deflate(fsv[i]);
		}
	}
	return ret;
}

void Report::strip_path_sfx(std::string& path, std::string sfx)
{
	auto fnb = std::filesystem::path(path).filename().string();
	auto dir = std::filesystem::path(path).parent_path();
	auto pos = fnb.find(sfx);
	if (pos != std::string::npos) {
		fnb.replace(pos, sfx.size(), "");
		auto fout = dir / fnb;
		INFO("moving ", path, " -> ", fout);
		std::filesystem::rename(path, fout);
	}
	else {
		WARN("no ", sfx, " found in ", path);
	}
}