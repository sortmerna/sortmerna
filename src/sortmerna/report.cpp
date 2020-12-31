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
	for (auto i = 0; i < vzlib_out.size(); ++i) {
		vzlib_out[i].init(true);
	}

	// prepare Readstates OUT
	vstate_out.resize(fv.size());
}

void Report::merge(int num_splits)
{
	// merge if there is more than one split
	if (fv.size() > 1) {
		std::ofstream ofs(fv[0], std::ios_base::app | std::ios_base::binary);
		if (!ofs.is_open()) {
			ERR("failed to open for writing: ", fv[0]);
			exit(1);
		}
		for (int i = 1; i < num_splits; ++i) {
			std::ifstream ifs(fv[i], std::ios_base::out | std::ios_base::binary);
			if (ifs.is_open()) {
				ofs << ifs.rdbuf();
				INFO("merged ", fv[i], " -> ", fv[0]);
				ifs.close();
				std::filesystem::remove(fv[i]);
				INFO("deleted ", fv[i]);
			}
			else {
				ERR("failed to open for reading: ", fv[i]);
				exit(1);
			}
		}
	}

	strip_path_sfx(fv[0]);
}

void Report::openfw(int idx)
{
	if (!fsv[idx].is_open()) {
		fsv[idx].open(fv[idx], std::ios::binary | std::ios::app);
	}
	if (!fsv[idx].good()) {
		ERR("Could not open output file [", fv[idx], "] for writing.");
		exit(EXIT_FAILURE);
	}
	else {
		INFO("Opened output file ", fv[idx], " for writing.");
	}
}

void Report::openfw()
{
	for (size_t i = 0; i < fv.size(); ++i) {
		openfw(i);
	}
}

void Report::openfr(int idx)
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

void Report::closef(int idx)
{
	if (fsv[idx].is_open()) {
		fsv[idx].flush();
		fsv[idx].close();
		INFO("Closed output file: ", fv[idx]);
	}
}

void Report::closef()
{
	for (int i = 0; i < fsv.size(); ++i) {
		closef(i);
	}
}

int Report::finish_deflate()
{
	int ret = 0;
	if (is_zip) {
		for (int i = 0; i < vzlib_out.size(); ++i) {
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