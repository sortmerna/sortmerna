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
#include <array>

#include "common.hpp"
#include "kvdb.hpp"

// global constants
const std::string \
OPT_REF = "ref",
OPT_READS = "reads",
OPT_ALIGNED = "aligned",
OPT_OTHER = "other",
OPT_WORKDIR = "workdir",
OPT_FASTX = "fastx",
OPT_SAM = "sam",
OPT_SQ = "SQ",
OPT_BLAST = "blast",
OPT_LOG = "log",
OPT_NUM_ALIGNMENTS = "num_alignments",
OPT_BEST = "best",
OPT_MIN_LIS = "min_lis",
OPT_PRINT_ALL_READS = "print_all_reads",
OPT_PAIRED_IN = "paired_in",
OPT_PAIRED_OUT = "paired_out",
OPT_MATCH = "match",
OPT_MISMATCH = "mismatch",
OPT_GAP_OPEN = "gap_open",
OPT_GAP_EXT = "gap_ext",
OPT_A = "a",
OPT_D = "d",
OPT_E = "e",
OPT_F = "F",
OPT_H = "h",
OPT_L = "L",
OPT_M = "m",
OPT_N = "N",
OPT_R = "R",
OPT_V = "v",
OPT_ID = "id",
OPT_COVERAGE = "coverage",
OPT_DE_NOVO_OTU = "de_novo_otu",
OPT_OTU_MAP = "otu_map",
OPT_PASSES = "passes",
OPT_EDGES = "edges",
OPT_NUM_SEEDS = "num_seeds",
OPT_FULL_SEARCH = "full_search",
OPT_PID = "pid",
OPT_VERSION = "version",
OPT_CMD = "cmd",
OPT_TASK = "task",
OPT_THREADS = "threads",
OPT_THPP = "thpp",
OPT_THREP = "threp",
OPT_DBG_PUT_DB = "dbg_put_db",
OPT_TMPDIR = "tmpdir",
OPT_INTERVAL = "interval",
OPT_MAX_POS = "max_pos";

// help strings
const std::string \
help_header =
"  Usage:   sortmerna --ref FILE [--ref FILE] --reads FWD_READS [--reads REV_READS] [OPTIONS]:\n"
"  -------------------------------------------------------------------------------------------------------------\n"
"  | option            type-format           description                                            default    |\n"
"  -------------------------------------------------------------------------------------------------------------\n",
help_ref = 
	"Reference file (FASTA) absolute or relative path.\n"
	"                                            Use mutliple times, once per a reference file",
help_reads = 
	"Raw reads file (FASTA/FASTQ).\n"
	"                                            Use twice for files with paired reads",
help_aligned = 
	"Aligned reads file name prefix.                         'aligned'\n"
	"                                            TODO: Remove",
help_other = 
	"Non-aligned reads output file name prefix               'other'\n"
	"                                            TODO: Remove.",
help_fastx = 
	"Output aligned reads into FASTA/FASTQ file",
help_workdir = 
	"Working directory path for storing Reference   USRDIR/sortmerna/\n"
	"                                            index, Key-value database and the output.\n",
help_sam = 
	"Output SAM alignment for aligned reads.",
help_SQ = 
	"Add SQ tags to the SAM file",
help_blast = 
	"output alignments in various Blast - like formats\n"
	"                                            '0'                    - pairwise\n"
	"                                            '1'                    - tabular(Blast - m 8 format)\n"
	"                                            '1 cigar'              - tabular + column for CIGAR\n"
	"                                            '1 cigar qcov'         - tabular + columns for CIGAR\n"
	"                                                                     and query coverage\n"
	"                                            '1 cigar qcov qstrand' - tabular + columns for CIGAR,\n"
	"                                                                     query coverage and strand",
help_dbg_put_db = 
	"",
help_log = 
	"Output overall statistics.                              True\n"
	"                                            TODO: remove",
help_num_alignments = 
	"Positive integer (INT >=0). Optional.\n"
	"                                            Report first INT alignments per read reaching E-value\n"
	"                                            If INT = 0, all alignments will be output",
help_best = 
	"Report INT best alignments per read reaching E-value    1\n"
	"                                            by searching --min_lis INT candidate alignments\n"
	"                                            (if 0 - all candidate alignments will be searched)",
help_min_lis = 
	"Search all alignments having the first INT longest LIS\n"
	"                                            LIS stands for Longest Increasing Subsequence,\n"
	"                                            it is computed using seeds' positions to expand hits into\n"
	"                                            longer matches prior to Smith - Waterman alignment.",
help_print_all_reads = 
	"Output null alignment strings for non-aligned reads     False\n"
	"                                            to SAM and/or BLAST tabular files",
help_paired_in = 
	"If one of the paired-end reads is Aligned,              False\n"
	"                                            put both reads into Aligned FASTA/Q file",
help_paired_out = 
	"If one of the paired-end reads is Non-aligned,          False\n"
	"                                            put both reads into Non-Aligned FASTA/Q file",
help_match = 
	"SW score (positive integer) for a match.                2",
help_mismatch = 
	"SW penalty (negative integer) for a mismatch.          -3",
help_gap_open = 
	"SW penalty (positive integer) for introducing a gap.    5",
help_gap_ext = 
	"SW penalty (positive integer) for extending a gap.      2",
help_N = 
	"SW penalty for ambiguous letters (N's) scored\n"
	"                                            as --mismatch",
help_F = 
	"Search only the forward strand.                         False",
help_R = 
	"Search only the reverse-complementary strand.           False",
help_e = 
	"E-value threshold.                                      1",
help_v = 
	"Produce verbose output.                                 False",
help_id = 
	"%%id similarity threshold (the alignment                0.97\n"
	"                                            must still pass the E-value threshold).",
help_coverage = 
	"%%query coverage threshold (the alignment must          0.97\n"
	"                                            still pass the E-value threshold)",
help_de_novo_otu = 
	"FASTA/FASTQ file for reads matching database < %%id     False\n"
	"                                            (set using --id) and < %%cov (set using --coverage)\n"
	"                                            (alignment must still pass the E-value threshold).",
help_otu_map = 
	"Output OTU map (input to QIIME's make_otu_table.py).    False",
help_passes = 
	"Three intervals at which to place the seed on the read  L,L/2,3\n"
	"                                            (L is the seed length)",
help_edges = 
	"Number (or percent if INT followed by %% sign) of       4\n"
	"                                            nucleotides to add to each edge of the read\n"
	"                                            prior to SW local alignment",
help_num_seeds = 
	"Number of seeds matched before searching                2\n"
	"                                            for candidate LIS",
help_pid = 
	"Add pid to output file names.                           False",
help_full_search = 
	"Search for all 0-error and 1-error seed                 False\n"
	"                                            matches in the index rather than stopping\n"
	"                                            after finding a 0-error match (<1%% gain in\n"
	"                                            sensitivity with up four-fold decrease in speed)",
help_h = 
	"Print help information",
help_version = 
	"Print SortMeRNA version number",
help_cmd = 
	"Launch an interactive session (command prompt)          False",
help_task = 
	"Processing Task:                                        4\n"
	"                                            0 - align. Only perform alignment\n"
	"                                            1 - post-processing (log writing)\n"
	"                                            2 - generate reports\n"
	"                                            3 - align and post-process\n"
	"                                            4 - all",
help_d 
	= "key-value datastore FULL folder path.              WORKDIR/kvdb/",
help_a = 
	"Number of threads to use                                numCores",
help_threads = 
	"Number of Read:Write:Process threads to use             1:1:numCores",
help_thpp = 
	"Number of Post-Processing Read:Process threads to use   1:1",
help_threp = 
	"Number of Report Read:Process threads to use            1:1",
help_tmpdir = 
	"Indexing: directory for writing temporary files when\n"
	"                                            building the reference index",
help_interval = 
	"Indexing: Positive integer: index every Nth L-mer in    1\n"
	"                                            the reference database e.g. '--interval 2'.",
help_m = 
	"Indexing: the amount of memory (in Mbytes) for building 3072\n"
	"                                            the index.",
help_L = 
	"Indexing: seed length.                                  18",
help_max_pos = 
	"Indexing: maximum (integer) number of positions to store  1000\n"
	"                                            for each unique L-mer. If 0 all positions are stored."
;

/* 
 * 1. 'blastops' 
 *     Vector of strings to store result from option --blast STRING.
 *    + --blast '0': output pairwise alignments\n
 *    + --blast '1': output BLAST Tabular format with the fields:
 *		   queryId, subjectId, percIdentity, alnLength, mismatchCount,
 *		   gapOpenCount, queryStart, queryEnd, subjectStart, subjectEnd, eVal, bitScore\n
 *    + --blast '1 cigar': tabular format + CIGAR string\n
 *    + --blast '1 cigar qcov': tabular format + CIGAR string + % query coverage\n
 *    + --blast '1 cigar qcov strand': tabular format + CIGAR string + % query coverage + strand\n
 * 2. 'skiplengths'
 *      '--passes' - for each index file three intervals at which to place the k-mer window on the read when searching for matches. 
 *      Defaults are calculated in Refstats::load e.g. {18,9,3} as follows:
 *
 *      ------------------------------------------------------	Read len = 54 (example). Max 13 positions to test
 *      1                 2                 3					Pass 1, step: 18 ------------------
 *      ^        4        ^        5        ^        			Pass 2, step: 9  ---------
 *      ^  6  7  ^ 	8  9  ^  |  |  ^  |  |  ^					Pass 3, step: 3  ---
 *
 * 3. 'num_alignments'
 *      unlike '--best', which searches many alignments(specified by '--min_lis') prior to outputting the best ones.
 * 4. 'align_cov'
 *      query coverage threshold (the alignment must still pass the E-value threshold)
 */
struct Runopts 
{
public:
	Runopts(int argc, char**argv, bool dryrun);
	~Runopts() {}

	enum OPT_CATEGORY { COMMON, OTU_PICKING, ADVANCED, DEVELOPER, HELP, INDEXING };

	typedef void (Runopts::*OptsMemFunc)(const std::string&); // pointer to member function
	typedef std::tuple<std::string, std::string, OPT_CATEGORY, bool, std::string, OptsMemFunc> opt_6_tuple;
	//                         |          |         |           |           |          |_ pointer to option processing function
	//                         |          |         |           |           |_Help string
	//                         |          |         |           |_Required option flag
	//                         |          |         |_Category of option
	//                         |          |_Type of option value
	//                         |_Name of option

	const std::map<OPT_CATEGORY, std::string> opt_category_name_map {
		{COMMON, "COMMON"},
		{OTU_PICKING, "OTU_PICKING"},
		{ADVANCED, "ADVANCED"},
		{DEVELOPER, "DEVELOPER"},
		{HELP, "HELP"},
		{INDEXING, "INDEXING"}
	};

	void print_help();

public:
	std::string workdir; // Directory for index, KVDB, Output
	std::string cmdline;
	std::string kvdbPath; // help_d TODO: remove
	std::string aligned_out_pfx = "aligned"; // '--aligned' aligned reads output file prefix
	std::string other_out_pfx = "other"; // '--other' rejected reads output file prefix

	int num_read_thread = 1; // number of threads reading the Reads file.
	int num_write_thread = 1; // number of threads writing to Key-value database
	int num_proc_thread = 0; // '-a' number of threads to use for alignment, post-processing, reporting. Default - all available cores.
	int num_read_thread_pp = 1; // number of post-processing read threads
	int num_proc_thread_pp = 1; // number of post-processing processor threads
	int num_read_thread_rep = 1; // number of report reader threads
	int num_proc_thread_rep = 1; // number of report processor threads

	int queue_size_max = 100; // max number of Reads in the Read and Write queues. 10 works OK.

	int32_t num_alignments = -1; // [3] help_num_alignments
	int32_t min_lis = -1; // '--min_lis N' search all alignments having the first N longest LIS
	int32_t seed_hits = -1;
	int32_t num_best_hits = 0;
	int32_t edges = -1;

	uint32_t minoccur = 0; // TODO: add to cmd options. Min number of k-mer occurrences in the DB to use for matching. See 'index.lookup_tbl[kmer_idx].count'

	long match = 2; // '--match' SW score (positive integer) for a match               TODO: change to int8_t
	long mismatch = -3; // '--mismatch' SW penalty (negative integer) for a mismatch   TODO: change to int8_t
	long gap_open = 5; // '--gap_open' SW penalty (positive integer) for introducing a gap
	long gap_extension = 2; // '--gap_ext' SW penalty (positive integer) for extending a gap
	long score_N = 0; // '-N' SW penalty for ambiguous letters (N's)                   TODO: change to int8_t

	double evalue = -1.0; // '-e' E-value threshold
	double align_id = -1.0; // OTU-picking option: minimum %%id to keep alignment
	double align_cov = -1.0; // [4] '--coverage': OTU-picking option: minimum %%coverage to keep alignment.

	// indexing options
	double max_file_size = 3072; // max size of an index file (or a part of the file). When exceeded, the index is split into parts.
	uint32_t lnwin_gv = 18;
	uint32_t interval = 1; // size of k-mer window shift. Default 1 is the min possible to generate max number of k-mers.
	uint32_t max_pos = 10000;
	// ~ END indexing options

	bool forward = false; // '-F' search only the forward strand if true
	bool reverse = false; // '-R' search only the reverse-complementary strand if true
	bool pairedin = false; // '--paired_in' both paired-end reads go in 'aligned' fasta/q file. Only Fasta/q and De-novo reporting.
	bool pairedout = false; // '--paired_out' both paired-end reads go in 'other' fasta/q file. Only Fasta/q and De-novo reporting.
	bool de_novo_otu = false; // '--de_novo_otu' FASTA/FASTQ file for reads matching database < %%id (set using --id) and < %%cov (set using --coverage)
	bool write_log = true; // '--log' output overall statistics. TODO: remove this option, always generate the log.
	bool print_all_reads = false; // '--print_all_reads' output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files
	bool samout = false; // '--sam' output SAM alignment (for aligned reads only)
	bool blastout = false; // '--blast' output alignments in various Blast-like formats
	bool is_fastxout = false; // '--fastx' output FASTA/FASTQ file (for aligned and/or rejected reads)
	bool otumapout = false; // '--otu_map' output OTU map (input to QIIME's make_otu_table.py)
	bool pid = false; // --pid add pid to output file names
	bool as_percent = false;
	bool full_search = false;
	bool exit_early = false; // flag to exit processing when either the reads or the reference file is empty or not FASTA/FASTQ
	bool verbose; // Verbose mode
	bool is_gz = false; // flags reads file is compressed and can be read
	bool yes_SQ = false; // --SQ add SQ tags to the SAM file
	bool interactive = false; // start interactive session
	bool is_index_built = false; // flags the index is built and ready for use
	bool is_other = false; // flags to produce 'other' files

	bool dbg_put_kvdb = false; // DEBUG option. if True - do Not put records into Key-value DB. Debugging Memory Consumption.

	std::vector<std::string> blastops; // [1]
	std::vector<std::string> readfiles; // '--reads'
	std::vector<std::pair<std::string, std::string>> indexfiles; // "--refs" Pairs (Reference file, Index name)
	std::vector<std::vector<uint32_t>> skiplengths; // [2] '--passes' K-mer window shift sizes. Refstats::load

public:
	std::string dbkey = "run_options";
	const std::string IDX_DIR  = "idx";
	const std::string KVDB_DIR = "kvdb";
	const std::string OUT_DIR  = "out";

	enum ALIGN_REPORT { align, postproc, report, alipost, all };
	ALIGN_REPORT alirep = all;
	BlastFormat blastFormat = BlastFormat::TABULAR;

private:
	// methods
	void process(int argc, char**argv, bool dryrun);
	void validate();
	void opt_sort();

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
	std::string to_string();
	std::string to_bin_string();
	void store_to_db(KeyValueDatabase &kvdb);

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
	bool have_reads = false; // flags reads file is plain text and can be read

	// container for options passed to the program
	std::multimap<std::string, std::string> mopt;

	// OPTIONS Map - specifies all possible options
	//std::map<std::string, std::tuple<bool, std::string, void(*)(const std::string&)>> options
	//
	//std::make_tuple(OPT_ALIGNED,        "STRING",      COMMON,      false, help_aligned, &Runopts::opt_aligned),
	//std::make_tuple(OPT_OTHER,          "STRING",      COMMON,      false, help_other, &Runopts::opt_other),
	//std::make_tuple(OPT_LOG,            "BOOL",        COMMON,      false, help_log, &Runopts::opt_log),
	//std::make_tuple(OPT_D,              "BOOL",        COMMON,      false, help_d, &Runopts::opt_d),
	//std::make_tuple(OPT_TMPDIR,         "INT",         INDEXING,    false, help_tmpdir, &Runopts::opt_tmpdir),
	//
	const std::array<opt_6_tuple, 44> options = {
		std::make_tuple(OPT_REF,            "PATH",        COMMON,      true,  help_ref, &Runopts::opt_ref),
		std::make_tuple(OPT_READS,          "PATH",        COMMON,      true,  help_reads, &Runopts::opt_reads),
		std::make_tuple(OPT_WORKDIR,        "PATH",        COMMON,      false, help_workdir, &Runopts::opt_workdir),
		std::make_tuple(OPT_FASTX,          "BOOL",        COMMON,      false, help_fastx, &Runopts::opt_fastx),
		std::make_tuple(OPT_SAM,            "BOOL",        COMMON,      false, help_sam, &Runopts::opt_sam),
		std::make_tuple(OPT_SQ,             "BOOL",        COMMON,      false, help_SQ, &Runopts::opt_SQ),
		std::make_tuple(OPT_BLAST,          "BOOL",        COMMON,      false, help_blast, &Runopts::opt_blast),
		std::make_tuple(OPT_NUM_ALIGNMENTS, "INT",         COMMON,      false, help_num_alignments, &Runopts::opt_num_alignments),
		std::make_tuple(OPT_BEST,           "INT",         COMMON,      false, help_best, &Runopts::opt_best),
		std::make_tuple(OPT_MIN_LIS,        "INT",         COMMON,      false, help_min_lis, &Runopts::opt_min_lis),
		std::make_tuple(OPT_PRINT_ALL_READS,"BOOL",        COMMON,      false, help_print_all_reads, &Runopts::opt_print_all_reads),
		std::make_tuple(OPT_PAIRED_IN,      "BOOL",        COMMON,      false, help_paired_in, &Runopts::opt_paired_in),
		std::make_tuple(OPT_PAIRED_OUT,     "BOOL",        COMMON,      false, help_paired_out, &Runopts::opt_paired_out),
		std::make_tuple(OPT_MATCH,          "INT",         COMMON,      false, help_match, &Runopts::opt_match),
		std::make_tuple(OPT_MISMATCH,       "INT",         COMMON,      false, help_mismatch, &Runopts::opt_mismatch),
		std::make_tuple(OPT_GAP_OPEN,       "INT",         COMMON,      false, help_gap_open, &Runopts::opt_gap_open),
		std::make_tuple(OPT_GAP_EXT,        "INT",         COMMON,      false, help_gap_ext, &Runopts::opt_gap_ext),
		std::make_tuple(OPT_A,              "INT",         COMMON,      false, help_a, &Runopts::opt_a),
		std::make_tuple(OPT_E,              "DOUBLE",      COMMON,      false, help_e, &Runopts::opt_e),
		std::make_tuple(OPT_F,              "BOOL",        COMMON,      false, help_F, &Runopts::opt_F),
		std::make_tuple(OPT_N,              "BOOL",        COMMON,      false, help_N, &Runopts::opt_N),
		std::make_tuple(OPT_R,              "BOOL",        COMMON,      false, help_R, &Runopts::opt_R),
		std::make_tuple(OPT_ID,             "INT",         OTU_PICKING, false, help_id, &Runopts::opt_id),
		std::make_tuple(OPT_COVERAGE,       "INT",         OTU_PICKING, false, help_coverage, &Runopts::opt_coverage),
		std::make_tuple(OPT_DE_NOVO_OTU,    "BOOL",        OTU_PICKING, false, help_de_novo_otu, &Runopts::opt_de_novo_otu),
		std::make_tuple(OPT_OTU_MAP,        "BOOL",        OTU_PICKING, false, help_otu_map, &Runopts::opt_otu_map),
		std::make_tuple(OPT_PASSES,         "BOOL",        ADVANCED,    false, help_passes, &Runopts::opt_passes),
		std::make_tuple(OPT_EDGES,          "BOOL",        ADVANCED,    false, help_edges, &Runopts::opt_edges),
		std::make_tuple(OPT_NUM_SEEDS,      "BOOL",        ADVANCED,    false, help_num_seeds, &Runopts::opt_num_seeds),
		std::make_tuple(OPT_FULL_SEARCH,    "INT",         ADVANCED,    false, help_full_search, &Runopts::opt_full_search),
		std::make_tuple(OPT_PID,            "BOOL",        ADVANCED,    false, help_pid, &Runopts::opt_pid),
		std::make_tuple(OPT_L,              "DOUBLE",      INDEXING,    false, help_L, &Runopts::opt_L),
		std::make_tuple(OPT_M,              "DOUBLE",      INDEXING,    false, help_m, &Runopts::opt_m),
		std::make_tuple(OPT_V,              "INT,INT,INT", INDEXING,    false, help_v, &Runopts::opt_v),
		std::make_tuple(OPT_INTERVAL,       "INT",         INDEXING,    false, help_interval, &Runopts::opt_interval),
		std::make_tuple(OPT_MAX_POS,        "INT",         INDEXING,    false, help_max_pos, &Runopts::opt_max_pos),
		std::make_tuple(OPT_H,              "BOOL",        HELP,        false, help_h, &Runopts::opt_h),
		std::make_tuple(OPT_VERSION,        "INT",         HELP,        false, help_version, &Runopts::opt_version),
		std::make_tuple(OPT_DBG_PUT_DB,     "BOOL",        DEVELOPER,   false, help_dbg_put_db, &Runopts::opt_dbg_put_db),
		std::make_tuple(OPT_CMD,            "INT:INT:INT", DEVELOPER,   false, help_cmd, &Runopts::opt_cmd),
		std::make_tuple(OPT_TASK,           "INT:INT:INT", DEVELOPER,   false, help_task, &Runopts::opt_task),
		std::make_tuple(OPT_THREADS,        "INT:INT:INT", DEVELOPER,   false, help_threads, &Runopts::opt_threads),
		std::make_tuple(OPT_THPP,           "BOOL",        DEVELOPER,   false, help_thpp, &Runopts::opt_thpp),
		std::make_tuple(OPT_THREP,          "PATH",        DEVELOPER,   false, help_threp, &Runopts::opt_threp)
	};
	// ~map options
}; // ~struct Runopts
// ~options.cpp
