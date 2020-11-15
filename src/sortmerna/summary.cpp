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
	pid_str = std::to_string(getpid());
	if (opts.is_pid)
	{
		sfx += "_" + pid_str;
	}
	sfx += ".log";
	f_log = opts.aligned_pfx.string() + sfx;
	INFO("Using summary file: ", f_log.generic_string());
	ofs_log.open(f_log, std::ofstream::binary | std::ofstream::out);
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

/*
 * Called for each index*index_part*read
 *
 * populate 'readstats.otu_map'
 * count 'readstats.total_reads_denovo_clustering'
 */
void writeSummary3(Read& read, Readstats& readstats, Refstats& refstats, References& refs, Runopts& opts)
{
	//uint32_t index_max_score = read.alignment.max_index; // index of alignment holding maximum SW score

	if (opts.is_otu_map || opts.is_denovo_otu) {
		if (read.is03) read.flip34();
	}

	// populate OTU map
	if (opts.is_otu_map) {
		if (read.is_id && read.is_cov) {
			// reference sequence identifier for mapped read
			std::string refhead = refs.buffer[read.alignment.alignv[read.alignment.max_index].ref_num].header;
			std::string ref_seq_str = refhead.substr(0, refhead.find(' '));
			// left trim '>' or '@'
			ref_seq_str.erase(ref_seq_str.begin(),
				std::find_if(ref_seq_str.begin(), ref_seq_str.end(),
					[](auto ch) {return !(ch == FASTA_HEADER_START || ch == FASTQ_HEADER_START);}));

			// read identifier
			std::string read_seq_str = read.getSeqId();
			readstats.pushOtuMap(ref_seq_str, read_seq_str); // thread safe
		}
	}

	// only call once per read, on the last index/part
	if (opts.is_denovo_otu) {
		if (refs.num == opts.indexfiles.size() - 1
			&& refs.part == refstats.num_index_parts[opts.indexfiles.size() - 1] - 1
			&& read.is_hit
			&& read.is_denovo)
		{
			++readstats.total_reads_denovo_clustering;
		}
	}
} // ~writeSummary3

/*
  runs in a thread
*/
void writeSummary2(int id, Readfeed& readfeed, Runopts& opts, References& refs, Readstats& readstats, Refstats& refstats, KeyValueDatabase& kvdb)
{
	unsigned countReads = 0;
	unsigned count_reads_aligned = 0;
	std::string readstr;

	INFO("PostProcessor: ", id, " thread: ", std::this_thread::get_id(), " started");

	for (;readfeed.next(id, readstr);)
	{
		{
			Read read(readstr);
			read.init(opts);
			read.load_db(kvdb);

			if (!read.isValid)
				continue;

			writeSummary3(read, readstats, refstats, refs, opts);
			readstr.resize(0);
			++countReads;
			if (read.is_hit) ++count_reads_aligned;
		}
	}

	INFO("PostProcessor: ", id, " thread: ", std::this_thread::get_id(),
		" done. Processed reads: ", countReads, ". count_reads_aligned: ", count_reads_aligned);
}

// called from main
void writeSummary(Readfeed& readfeed, Readstats& readstats, KeyValueDatabase& kvdb, Runopts& opts)
{
	INFO("==== Starting summary (alignment statistics report) ====\n\n");
	Summary summary;
	bool indb = readstats.restoreFromDb(kvdb);
	if (indb)
		INFO("Restored Readstats from DB:\n    ", readstats.toString());

	readstats.total_reads_denovo_clustering = 0; // TODO: to prevent incrementing the stored value. Change this if ever using 'stats_calc_done"

	//if (!readstats.stats_calc_done)
	//{
	Refstats refstats(opts, readstats);

	// this part is only necessary for OTU map and/or deNovo clustering
	if (opts.is_otu_map || opts.is_denovo_otu) {
		//ReadsQueue read_queue("queue_1", opts.queue_size_max, readstats.all_reads_count);
		int numThreads = 0;
		if (opts.feed_type == FEED_TYPE::LOCKLESS)
		{
			numThreads = opts.num_read_thread_pp + opts.num_proc_thread_pp;
			INFO("using total threads: ", numThreads, " including Read threads: ", opts.num_read_thread_pp, " Processor threads: ", opts.num_proc_thread_pp);
		}
		else {
			numThreads = opts.num_proc_thread_pp;
			INFO("Using total threads: ", numThreads);
		}

		std::vector<std::thread> tpool;
		tpool.reserve(numThreads);

		References refs;
		// loop through every reference file part
		for (uint16_t idx = 0; idx < opts.indexfiles.size(); ++idx) {
			for (uint16_t ipart = 0; ipart < refstats.num_index_parts[idx]; ++ipart) {
				// load reference
				INFO("Loading reference ", idx, " part ", ipart + 1, "/", refstats.num_index_parts[idx], "  ... ");
				auto starts = std::chrono::high_resolution_clock::now();
				refs.load(idx, ipart, opts, refstats);
				std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - starts;
				INFO("done. Elapsed sec: [", elapsed.count(), "]");

				starts = std::chrono::high_resolution_clock::now(); // index processing starts

				// add Readfeed job if necessary
				if (opts.feed_type == FEED_TYPE::LOCKLESS)
				{
					//tpool.addJob(f_readfeed_run);
					//tpool.addJob(Readfeed(opts.feed_type, opts.readfiles, opts.is_gz));
				}
				tpool.emplace_back(std::thread(writeSummary2, 0, std::ref(readfeed), std::ref(opts), std::ref(refs),
					std::ref(readstats), std::ref(refstats), std::ref(kvdb)));

				// wait till processing is done on one index part
				//tpool.waitAll(); 
				for (auto i = 0; i < tpool.size(); ++i) {
					tpool[i].join();
				}

				refs.unload();
				//read_queue.reset();

				elapsed = std::chrono::high_resolution_clock::now() - starts;
				INFO_MEM("Done reference ", idx, " Part: ", ipart + 1, " Elapsed sec: ", elapsed.count());
			} // ~for(ipart)
		} // ~for(idx)

		INFO("total_reads_denovo_clustering = ", readstats.total_reads_denovo_clustering);
	} // ~if opts.is_otu_map || opts.is_denovo_otu

	readstats.set_is_total_mapped_sw_id_cov();
	readstats.is_stats_calc = true;
	readstats.store_to_db(kvdb); // store reads statistics computed by post-processor
//} // ~if !readstats.stats_calc_done

	summary.write(opts, refstats, readstats);

	if (opts.is_otu_map) {
		// OTU map output file  WORKDIR/out/aligned_otus.txt
		std::string f_otumap;
		std::ofstream otumap;
		std::string sfx;
		if (opts.is_pid)
		{
			sfx += "_" + std::to_string(getpid());
		}
		sfx += "_otus.txt";
		f_otumap = opts.aligned_pfx.string() + sfx;
		readstats.printOtuMap(f_otumap);
	}

	INFO("==== Done Post-processing (alignment statistics report) ====\n\n");
} // ~writeSummary