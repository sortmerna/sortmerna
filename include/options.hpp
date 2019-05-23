#pragma once
/**
 * FILE: options.hpp
 * Created: Aug 19, 2017 Sat
 *
 * skiplength
 *    skip lengths for pass 1, pass 2 and pass 3 in first step of sortmerna
 *    pipeline for each reference database searched
 */
#include <string>
#include <vector>
#include <map>
#include <tuple>

#include "common.hpp"

struct Runopts 
{
	Runopts(int argc, char**argv, bool dryrun);
	~Runopts() {}

	typedef void (Runopts::*OptsMemFunc)(const std::string&); // pointer to member function

	std::string kvdbPath; // '-d' (opt_d) key-value database for alignment results
	std::vector<std::string> readfiles;
	std::string readsrev; // '--reads' | '--reads-gz' reversed reads file when processing paired reads
	std::string filetype_ar; // '--aligned' aligned reads output file
	std::string filetype_or; // '--other' rejected reads output file
	std::string cmdline;
	std::string workdir;
	const std::string IDX_DIR  = "idx";
	const std::string KVDB_DIR = "kvdb";
	const std::string OUT_DIR  = "out";

	int num_read_thread = 1; // number of threads reading the Reads file.
	int num_write_thread = 1; // number of threads writing to Key-value database
	int num_proc_thread = 0; // '-a' number of threads to use for alignment, post-processing, reporting. Default - all available cores.

	int num_read_thread_pp = 1; // number of post-processing read threads
	int num_proc_thread_pp = 1; // number of post-processing processor threads

	int num_read_thread_rep = 1; // number of report reader threads
	int num_proc_thread_rep = 1; // number of report processor threads

	int queue_size_max = 100; // max number of Reads in the Read and Write queues. 10 works OK.

	long match = 2; // '--match' SW score (positive integer) for a match               TODO: change to int8_t
	long mismatch = -3; // '--mismatch' SW penalty (negative integer) for a mismatch   TODO: change to int8_t
	long gap_open = 5; // '--gap_open' SW penalty (positive integer) for introducing a gap
	long gap_extension = 2; // '--gap_ext' SW penalty (positive integer) for extending a gap
	long score_N = 0; // '-N' SW penalty for ambiguous letters (N's)                   TODO: change to int8_t

	double evalue = -1.0; /* '-e' E-value threshold */
	double align_id = -1.0; /* OTU-picking option: minimum %%id to keep alignment */
	/* '--coverage' query coverage threshold (the alignment must still pass the E-value threshold) 
		OTU-picking option: minimum %%coverage to keep alignment. */
	double align_cov = -1.0;

	/* '--num_alignments': output the first '--num_alignments' found, unlike '--best', 
		which searches many alignments(specified by '--min_lis') prior to outputting the best ones. */
	int32_t num_alignments = -1; // 
	int32_t min_lis = -1; // '--min_lis' search all alignments having the first N longest LIS
	int32_t seed_hits = -1;
	int32_t num_best_hits = 0;
	int32_t edges = -1;

	uint32_t minoccur = 0; // TODO: add to cmd options. Min number of k-mer occurrences in the DB to use for matching. See 'index.lookup_tbl[kmer_idx].count'

	// indexing options
	double max_file_size = 3072; // max size of an index file (or a part of the file). When exceeded, the index is split into parts.
	uint32_t lnwin_gv = 18;
	uint32_t interval = 1;
	uint32_t max_pos = 10000;

	bool forward = false; // '-F' search only the forward strand if true
	bool reverse = false; // '-R' search only the reverse-complementary strand if true
	bool paired = false; // '--paired' flags processing paired reads
	bool pairedin = false; // '--paired_in' both paired-end reads go in --aligned fasta/q file. Only Fasta/q and De-novo reporting.
	bool pairedout = false; // '--paired_out' both paired-end reads go in --other fasta/q file. Only Fasta/q and De-novo reporting.
	bool de_novo_otu = false; // '--de_novo_otu' FASTA/FASTQ file for reads matching database < %%id (set using --id) and < %%cov (set using --coverage)
	bool doLog = false; // '--log' output overall statistics
	bool print_all_reads = false; // '--print_all_reads' output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files
	bool samout = false; // '--sam' output SAM alignment (for aligned reads only)
	bool blastout = false; // '--blast' output alignments in various Blast-like formats
	bool fastxout = false; // '--fastx' output FASTA/FASTQ file (for aligned and/or rejected reads)
	bool otumapout = false; // '--otu_map' output OTU map (input to QIIME's make_otu_table.py)
	bool pid = false;
	bool as_percent = false;
	bool full_search = false;
	bool exit_early = false; // flag to exit processing when either the reads or the reference file is empty or not FASTA/FASTQ
	bool verbose; // Verbose mode
	bool is_gz = false; // flags reads file is compressed and can be read
	bool yes_SQ = false; // --SQ add SQ tags to the SAM file
	bool interactive = false; // start interactive session
	bool is_index_built = false; // flags the index is built and ready for use

	// DEBUG options
	bool dbg_put_kvdb = false; // if True - do Not put records into Key-value DB. Debugging Memory Consumption.

	enum ALIGN_REPORT { align, postproc, report, alipost, all };
	ALIGN_REPORT alirep = all;
	BlastFormat blastFormat = BlastFormat::TABULAR;

	/*! @brief Vector of strings to store result from option --blast STRING.
		+ --blast '0': output pairwise alignments\n
		+ --blast '1': output BLAST Tabular format with the fields:
				queryId, subjectId, percIdentity, alnLength, mismatchCount, 
				gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore\n
		+ --blast '1 cigar': tabular format + CIGAR string\n
		+ --blast '1 cigar qcov': tabular format + CIGAR string + % query coverage\n
		+ --blast '1 cigar qcov strand': tabular format + CIGAR string + % query coverage + strand\n
	*/
	std::vector<std::string> blastops;
	// "--refs" Pairs (Reference file, Index name)
	std::vector<std::pair<std::string, std::string>> indexfiles;
	/* '--passes' (optional): for each index file three intervals at which to place the seed on the read. <-- Refstats::load
		Defaults: 0 */
	std::vector<std::vector<uint32_t>> skiplengths;

private:
	// methods
	void process(int argc, char**argv, bool dryrun);
	void validate();

	void opt_reads(const std::string &val);
	void opt_reads_gz(char **argv, int &narg);
	void opt_ref(const std::string &val);
	void opt_aligned(const std::string &val); // TODO: make optional. Aligned will be named automatically and put into WORKDIR
	void opt_other(const std::string &val); // TODO: make optional. Similar to Aligned.
	void opt_log(const std::string &val);
	void opt_de_novo_otu(const std::string &val);
	void opt_otu_map(const std::string &val);
	void opt_print_all_reads(const std::string &val);
	void opt_pid(const std::string &val);
	void opt_paired_in(const std::string &val);
	void opt_paired_out(const std::string &val);
	void opt_match(const std::string &val);
	void opt_mismatch(const std::string &val);
	void opt_gap_open(const std::string &val);
	void opt_gap_ext(const std::string &val);
	void opt_num_seeds(const std::string &val);
	void opt_fastx(const std::string &val);
	void opt_sam(const std::string &val);
	void opt_blast(const std::string &val);
	void opt_min_lis(const std::string &val);
	void opt_best(const std::string &val);
	void opt_num_alignments(const std::string &val);
	void opt_edges(const std::string &val);
	void opt_full_search(const std::string &val);
	void opt_SQ(const std::string &val);
	void opt_passes(const std::string &val);
	void opt_id(const std::string &val);
	void opt_coverage(const std::string &val);
	void opt_version(const std::string &val);
	void opt_task(const std::string &val);
	void opt_cmd(const std::string &val);
	void opt_threads(const std::string &val);
	void opt_thpp(const std::string &val); // post-proc threads --thpp 1:1
	void opt_threp(const std::string &val); // report threads --threp 1:1 
	void opt_a(const std::string &val);
	void opt_e(const std::string &val); // opt_e_Evalue
	void opt_F(const std::string &val); // opt_F_ForwardOnly
	void opt_R(const std::string &val); // opt_R_ReverseOnly
	void opt_h(const std::string &val);
	void opt_v(const std::string &val); // opt_v_Verbose
	void opt_N(const std::string &val); // opt_N_MatchAmbiguous
	void opt_d(const std::string &val); // opt_d_KeyValDatabase Key-Value Database directory path (kvdbPath)
	void opt_workdir(const std::string &path);

	// ref tmpdir interval m L max_pos v h  // indexing options
	void opt_tmpdir(const std::string &val);
	void opt_interval(const std::string &val);
	void opt_m(const std::string &val);
	void opt_L(const std::string &val);
	void opt_max_pos(const std::string &val);

	void opt_default(const std::string &opt);
	void opt_dbg_put_db(const std::string &opt);
	void opt_unknown(char **argv, int &narg, char * opt);

	void test_kvdb_path();

private:
	// SW alignment parameters
	bool match_set = false;
	bool mismatch_set = false;
	bool gap_open_set = false;
	bool gap_ext_set = false;
	bool full_search_set = false;
	bool passes_set = false;
	bool edges_set = false;
	bool match_ambiguous_N = false; // -N flags to match the ambiguous characters using score_N
	bool min_lis_set = false;
	bool num_alignments_set = false;
	bool best_set = false;
	bool have_reads = false; // '--reads' flags reads file is plain text and can be read

	// help strings
	std::string \
		help_ref = "FASTA reference file. Use mutliple '--ref' options to specify multiple files",
		help_reads = "FASTA/FASTQ raw reads file. Use multiple '--reads' options for multiple files",
		help_aligned = "Aligned reads file",
		help_other = "",
		help_fastx = "",
		help_workdir = "Working directory path where to put KVDB and all the output files.",
		help_sam = "",
		help_SQ = "",
		help_blast = "",
		help_dbg_put_db = "",
		help_log = "",
		help_num_alignments = "",
		help_best = "",
		help_min_lis = "",
		help_print_all_reads = "",
		help_paired_in = "",
		help_paired_out = "",
		help_match = "",
		help_mismatch = "",
		help_gap_open = "",
		help_gap_ext = "",
		help_N = "",
		help_F = "",
		help_R = "",
		help_e = "",
		help_v = "",
		help_id = "",
		help_coverage = "",
		help_de_novo_otu = "",
		help_otu_map = "",
		help_passes = "",
		help_edges = "",
		help_num_seeds = "",
		help_pid = "",
		help_full_search = "",
		help_h = "",
		help_version = "",
		help_cmd = "",
		help_task = "",
		help_d = "key-value datastore FULL folder path. Default: USERDIR/kvdb",
		help_a = "",
		help_threads = "",
		help_thpp = "",
		help_threp = "",
		help_tmpdir = "Indexing: directory for writing temporary files when building the reference index",
		help_interval = "Indexing: Positive integer: index every Nth L-mer in the reference database e.g. '--interval 2'. Default 1",
		help_m = "Indexing: the amount of memory (in Mbytes) for building the index. Default 3072",
		help_L = "Indexing: seed length. Default 18",
		help_max_pos = "Indexing: maximum (integer) number of positions to store for each unique L-mer. Default 1000. If 0 all positions are stored."
		;

	// container for options passed to the program
	std::multimap<std::string, std::string> mopt;

	// OPTIONS Map specifies all possible options
	//std::map<std::string, std::tuple<bool, std::string, void(*)(const std::string&)>> options
	std::map<std::string, std::tuple<bool, std::string, OptsMemFunc>> options
	{
		//     |                      |         |               |_pointer to option processing function
		//     |                      |         |_option help string
		//     |_option name          |_flag is option required
		{"ref",             {true, help_ref, &Runopts::opt_ref}},
		{"reads",           {true, help_reads, &Runopts::opt_reads}},
		{"aligned",         {false, help_aligned, &Runopts::opt_aligned}},
		{"other",           {false, help_other, &Runopts::opt_other}},
		{"workdir",         {false, help_workdir, &Runopts::opt_workdir}},
		{"fastx",           {false, help_fastx, &Runopts::opt_fastx}},
		{"sam",             {false, help_sam, &Runopts::opt_sam}},
		{"SQ",              {false, help_SQ, &Runopts::opt_SQ}},
		{"blast",           {false, help_blast, &Runopts::opt_blast}},
		{"log",             {false, help_log, &Runopts::opt_log}},
		{"num_alignments",  {false, help_num_alignments, &Runopts::opt_num_alignments}},
		{"best",            {false, help_best, &Runopts::opt_best}},
		{"min_lis",         {false, help_min_lis, &Runopts::opt_min_lis}},
		{"print_all_reads", {false, help_print_all_reads, &Runopts::opt_print_all_reads}},
		{"paired_in",       {false, help_paired_in, &Runopts::opt_paired_in}},
		{"paired_out",      {false, help_paired_out, &Runopts::opt_paired_out}},
		{"match",           {false, help_match, &Runopts::opt_match}},
		{"mismatch",        {false, help_mismatch, &Runopts::opt_mismatch}},
		{"gap_open",        {false, help_gap_open, &Runopts::opt_gap_open}},
		{"gap_ext",         {false, help_gap_ext, &Runopts::opt_gap_ext}},
		{"N",               {false, help_N, &Runopts::opt_N}},
		{"F",               {false, help_F, &Runopts::opt_F}},
		{"R",               {false, help_R, &Runopts::opt_R}},
		{"e",               {false, help_e, &Runopts::opt_e}},
		{"v",               {false, help_v, &Runopts::opt_v}},
		{"id",              {false, help_id, &Runopts::opt_id}},
		{"coverage",        {false, help_coverage, &Runopts::opt_coverage}},
		{"de_novo_otu",     {false, help_de_novo_otu, &Runopts::opt_de_novo_otu}},
		{"otu_map",         {false, help_otu_map, &Runopts::opt_otu_map}},
		{"passes",          {false, help_passes, &Runopts::opt_passes}},
		{"edges",           {false, help_edges, &Runopts::opt_edges}},
		{"num_seeds",       {false, help_num_seeds, &Runopts::opt_num_seeds}},
		{"full_search",     {false, help_full_search, &Runopts::opt_full_search}},
		{"pid",             {false, help_pid, &Runopts::opt_pid}},
		{"h",               {false, help_h, &Runopts::opt_h}},
		{"version",         {false, help_version, &Runopts::opt_version}},
		{"cmd",             {false, help_cmd, &Runopts::opt_cmd}},
		{"task",            {false, help_task, &Runopts::opt_task}},
		{"d",               {false, help_d, &Runopts::opt_d}},
		{"a",               {false, help_a, &Runopts::opt_a}},
		{"threads",         {false, help_threads, &Runopts::opt_threads}},
		{"thpp",            {false, help_thpp, &Runopts::opt_thpp}},
		{"threp",           {false, help_threp, &Runopts::opt_threp}},
		{"dbg_put_db",      {false, help_dbg_put_db, &Runopts::opt_dbg_put_db}},
		{"tmpdir",          {false, help_tmpdir, &Runopts::opt_tmpdir}},
		{"interval",        {false, help_interval, &Runopts::opt_interval}},
		{"m",               {false, help_m, &Runopts::opt_m}},
		{"L",               {false, help_L, &Runopts::opt_L}},
		{"max_pos",         {false, help_max_pos, &Runopts::opt_max_pos}}
	}; // ~map options
}; // ~struct Runopts
// ~options.cpp
