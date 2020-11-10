#include <fstream>

#include "summary.hpp"
#include "options.hpp"
#include "readstats.hpp"
#include "refstats.hpp"
#include "common.hpp"

Summary::Summary() :
	is_de_novo_otu(false),
	is_otumapout(false),
	total_reads(0),
	total_reads_denovo_clustering(0),
	total_reads_mapped(0),
	total_mapped_sw_id_cov(0),
	min_read_len(0),
	max_read_len(0),
	all_reads_len(0),
	total_otu(0)
{}

void Summary::write(Runopts& opts, Refstats& refstats, Readstats& readstats)
{
	std::filesystem::path f_log;
	std::ofstream ofs_log;
	std::string sfx;
	if (opts.is_pid)
	{
		sfx += "_" + std::to_string(getpid());
	}
	sfx += ".log";
	f_log = opts.aligned_pfx.string() + sfx;
	INFO("Using summary file: ", f_log.generic_string());
	ofs_log.open(f_log, std::ofstream::binary | std::ofstream::app);
	if (!ofs_log.is_open()) {
		ERR("Failed opening file ", f_log);
		exit(EXIT_FAILURE);
	}

	cmd = opts.cmdline;
	total_reads = readstats.all_reads_count;
	if (opts.is_denovo_otu) {
		is_de_novo_otu = opts.is_denovo_otu;
		total_reads_denovo_clustering = readstats.total_reads_denovo_clustering;
	}
	total_reads_mapped = readstats.total_reads_aligned.load(std::memory_order_relaxed);
	min_read_len = readstats.min_read_len;
	max_read_len = readstats.max_read_len;
	all_reads_len = readstats.all_reads_len;

	// stats by database
	for (uint32_t index_num = 0; index_num < opts.indexfiles.size(); index_num++)
	{
		auto pcn = (float)((float)readstats.reads_matched_per_db[index_num] / (float)readstats.all_reads_count) * 100;
		db_matches.emplace_back(std::make_pair(opts.indexfiles[index_num].first, pcn));
	}

	if (opts.is_otu_map) {
		is_otumapout = opts.is_otu_map;
		total_mapped_sw_id_cov = readstats.total_mapped_sw_id_cov.load(std::memory_order_relaxed);
		total_otu = readstats.otu_map.size();
	}

	// set timestamp  <ctime>
	std::time_t tm = std::time(0);
	timestamp = std::ctime(&tm); // Tue Oct 20 08:39:35 2020  Win deprecation: use 'ctime_s' or _CRT_SECURE_NO_WARNINGS

	ofs_log << to_string(opts, refstats);
	ofs_log.close();
} // ~Summary::write

std::string Summary::to_string(Runopts& opts, Refstats& refstats)
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
	if (is_de_novo_otu)
	{
		// all reads that have read::hit_denovo == true
		ss << "    Total reads for de novo clustering = " << total_reads_denovo_clustering << std::endl;
	}
	// output total non-rrna + rrna reads
	ss << std::setprecision(2) << std::fixed
		<< "    Total reads passing E-value threshold = " << total_reads_mapped
		<< " (" << ((float)total_reads_mapped / (float)total_reads * 100) << ")" << std::endl
		<< "    Total reads failing E-value threshold = " << total_reads - total_reads_mapped
		<< " (" << (1 - ((float)((float)total_reads_mapped / (float)total_reads))) * 100 << ")" << std::endl
		<< "    Minimum read length = " << min_read_len << std::endl
		<< "    Maximum read length = " << max_read_len << std::endl
		<< "    Mean read length    = " << all_reads_len / total_reads << std::endl << std::endl;

	ss << " Coverage by database:" << std::endl;

	// output stats by database
	for (auto match : db_matches)
	{
		ss << "    " << match.first << "\t\t" << match.second << std::endl;
	}

	if (is_otumapout)
	{
		ss << " Total reads passing %%id and %%coverage thresholds = " << total_mapped_sw_id_cov << std::endl
			<< " Total OTUs = " << total_otu << std::endl;
	}

	ss << std::endl << " " << timestamp << std::endl;

	return ss.str();
} // ~Summary::to_string