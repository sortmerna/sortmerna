#include <fstream>
#include <thread>

#include "summary.hpp"
#include "options.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "common.hpp"
#include "read.hpp"
#include "references.hpp"
#include "readfeed.hpp"

Summary::Summary() :
	is_de_novo(false),
	is_otumapout(false),
	total_reads(0),
	total_mapped(0),
	total_denovo(0),
	total_id_cov(0),
	total_otu(0),
	min_read_len(0),
	max_read_len(0),
	all_reads_len(0)
{}

void Summary::write(Refstats& refstats, Readstats& readstats, Runopts& opts)
{
	std::ofstream ofs;
	std::filesystem::path f_log;
	std::string sfx = opts.is_pid ? "_" + std::to_string(getpid()) : "";
	f_log = opts.aligned_pfx.string() + sfx + ".log";
	INFO("Using summary file: ", f_log.generic_string());
	ofs.open(f_log, std::ofstream::binary | std::ofstream::out);
	if (!ofs.is_open()) {
		ERR("Failed opening file ", f_log);
		exit(EXIT_FAILURE);
	}

	cmd = opts.cmdline;
	total_reads = readstats.all_reads_count;
	if (opts.is_denovo) {
		is_de_novo = opts.is_denovo;
		total_denovo = readstats.total_denovo;
	}
	total_mapped = readstats.total_aligned.load(std::memory_order_relaxed);
	min_read_len = readstats.min_read_len;
	max_read_len = readstats.max_read_len;
	all_reads_len = readstats.all_reads_len;

	// stats by database
	for (uint32_t index_num = 0; index_num < opts.indexfiles.size(); index_num++) {
		auto pcn = (float)((float)readstats.reads_matched_per_db[index_num] / (float)readstats.all_reads_count) * 100;
		db_matches.emplace_back(std::make_pair(opts.indexfiles[index_num].first, pcn));
	}

	if (opts.is_otu_map) {
		is_otumapout = opts.is_otu_map;
		total_id_cov = readstats.total_aligned_id_cov.load(std::memory_order_relaxed);
		total_otu = readstats.total_otu;
	}

	// set timestamp  <ctime>
	std::time_t tm = std::time(0);
	timestamp = std::ctime(&tm); // Tue Oct 20 08:39:35 2020  Win deprecation: use 'ctime_s' or _CRT_SECURE_NO_WARNINGS

	ofs << to_string(refstats, opts);
	ofs.close();
} // ~Summary::write

std::string Summary::to_string(Refstats& refstats, Runopts& opts)
{
	std::stringstream ss;
	size_t idx = 0;

	ss << " Command:\n    " << cmd << std::endl << std::endl

		<< " Process pid = " << pid_str << std::endl << std::endl

		<< " Parameters summary: " << std::endl;

	for (auto ref : opts.indexfiles) {
		ss << "    Reference file: " << ref.first << std::endl
			<< "        Seed length = " << opts.seed_win_len << std::endl
			<< "        Pass 1 = " << opts.skiplengths[idx][0]
			<< ", Pass 2 = " << opts.skiplengths[idx][1]
			<< ", Pass 3 = " << opts.skiplengths[idx][2] << std::endl
			<< "        Gumbel lambda = " << refstats.gumbel[idx].first << std::endl
			<< "        Gumbel K = " << refstats.gumbel[idx].second << std::endl
			<< "        Minimal SW score based on E-value = " << refstats.minimal_score[idx] << std::endl;
		++idx;
	}
	ss << "    Number of seeds = " << opts.hit_seeds << std::endl
		<< "    Edges = " << opts.edges << std::endl
		<< "    SW match = " << opts.match << std::endl
		<< "    SW mismatch = " << opts.mismatch << std::endl
		<< "    SW gap open penalty = " << opts.gap_open << std::endl
		<< "    SW gap extend penalty = " << opts.gap_extension << std::endl
		<< "    SW ambiguous nucleotide = " << opts.score_N << std::endl
		<< "    SQ tags are " << (opts.is_SQ ? "" : "not ") << "output" << std::endl
		<< "    Number of alignment processing threads = " << opts.num_proc_thread << std::endl;
	for (auto readf : opts.readfiles) {
		ss << "    Reads file: " << readf << std::endl;
	}
	ss << "    Total reads = " << total_reads << std::endl << std::endl;

	ss << " Results:" << std::endl;
	if (is_de_novo)
	{
		// all reads that have read::hit_denovo == true
		ss << "    Total reads for de novo clustering = " << total_denovo << std::endl;
	}
	// output total non-rrna + rrna reads
	auto ev_pass_ratio = (float)total_mapped / total_reads;
	ss << std::setprecision(2) << std::fixed
		<< "    Total reads passing E-value threshold = " << total_mapped
		<< " (" << (ev_pass_ratio * 100) << ")" << std::endl
		<< "    Total reads failing E-value threshold = " << total_reads - total_mapped
		<< " (" << (1 - ev_pass_ratio) * 100 << ")" << std::endl;

	if (is_otumapout)
	{
		auto idcov_pass_ratio = (float)total_id_cov / total_reads;
		ss << "    Total reads passing %%id and %%coverage thresholds = " << total_id_cov
			<< " (" << (idcov_pass_ratio * 100) << ")" << std::endl
		   << "    Total OTUs = " << total_otu << std::endl;
	}

	ss	<< "    Minimum read length = " << min_read_len << std::endl
		<< "    Maximum read length = " << max_read_len << std::endl
		<< "    Mean read length    = " << all_reads_len / total_reads << std::endl << std::endl;

	ss << " Coverage by database:" << std::endl;

	// output stats by database
	for (auto match : db_matches)
	{
		ss << "    " << match.first << "\t\t" << match.second << std::endl;
	}

	ss << std::endl << " " << timestamp << std::endl;

	return ss.str();
} // ~Summary::to_string


// called from main
void writeSummary(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	INFO("==== Starting summary of alignment statistics ====");
	auto start = std::chrono::high_resolution_clock::now();
	Refstats refstats(opts, readstats);
	Summary summary;
	summary.write(refstats, readstats, opts);
	std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
	INFO("==== Done summary in sec [", elapsed.count(), "] ====\n\n");
} // ~writeSummary