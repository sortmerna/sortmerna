#pragma once
/**
 * FILE: options.hpp
 * Created: Aug 19, 2017 Sat
 *
 * skiplength
 *    skip lengths for pass 1, pass 2 and pass 3 in first step of sortmerna
 *    pipeline for each reference database searched
 *
 * @copyright 2016-20 Clarity Genomics BVBA
 */
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <array>
#include <filesystem>

#include "common.hpp"
#include "kvdb.hpp"

// global constants
const std::string \
OPT_REF = "ref",
OPT_READS = "reads",
OPT_ALIGNED = "aligned",
OPT_OTHER = "other",
OPT_WORKDIR = "workdir",
OPT_KVDB = "kvdb",
OPT_IDX = "idx",
OPT_FASTX = "fastx",
OPT_SAM = "sam",
OPT_SQ = "SQ",
OPT_BLAST = "blast",
OPT_LOG = "log",
OPT_NUM_ALIGNMENTS = "num_alignments",
OPT_NO_BEST = "no-best",
OPT_MIN_LIS = "min_lis",
OPT_PRINT_ALL_READS = "print_all_reads",
OPT_PAIRED = "paired",
OPT_PAIRED_IN = "paired_in",
OPT_PAIRED_OUT = "paired_out",
OPT_OUT2 = "out2",
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
OPT_DENOVO_OTU = "de_novo_otu",
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
	"                                            Use mutliple times, once per a reference file\n",
help_reads = 
	"Raw reads file (FASTA/FASTQ/FASTA.GZ/FASTQ.GZ).\n"
	"                                            Use twice for files with paired reads\n",
help_aligned = 
	"Aligned reads file prefix [dir/][pfx]       WORKDIR/out/aligned\n"
	"                                            Directory and file prefix for aligned output i.e.\n"
	"                                            each output file will go into the specified directory with the given prefix.\n"
	"                                            The appropriate extension (fasta|fastq|blast|sam|etc) will be automatically added.\n"
	"                                            Both 'dir' and 'pfx' are optional.\n"
	"                                            The 'dir' can be a relative or an absolute path.\n"
	"                                            If 'dir' is not specified, the output will be created in the WORKDIR/out/\n"
	"                                            If 'pfx' is not specified, the prefix 'aligned' will be used\n"
	"                                            Examples:\n"
	"                                             -aligned $MYDIR/dir_1/dir_2/1 -> $MYDIR/dir_1/dir_2/1.fasta\n"
	"                                             -aligned dir_1/apfx           -> $PWD/dir_1/apfx.fasta\n"
	"                                             -aligned dir_1/               -> $PWD/aligned.fasta\n"
	"                                             -aligned apfx                 -> $PWD/apfx.fasta\n"
	"                                             -aligned  (no argument)       -> WORKDIR/out/aligned.fasta\n",
help_other = 
	"Non-aligned reads file prefix [dir/][pfx]   WORKDIR/out/other\n"
	"                                            Must be used with '" + OPT_FASTX + "'.\n"
	"                                            Directory and file prefix for non-aligned output i.e.\n"
	"                                            each output file will go into the specified directory with the given prefix.\n"
	"                                            The appropriate extension (fasta|fastq|blast|sam|etc) will be automatically added.\n"
	"                                            Both 'dir' and 'pfx' are optional.\n"
	"                                            The 'dir' can be a relative or an absolute path.\n"
	"                                            If 'dir' is not specified, the output will be created in the WORKDIR/out/\n"
	"                                            If 'pfx' is not specified, the prefix 'other' will be used\n"
	"                                            Examples:\n"
	"                                             -other $MYDIR/dir_1/dir_2/1 -> $MYDIR/dir_1/dir_2/1.fasta\n"
	"                                             -other dir_1/apfx           -> $PWD/dir_1/apfx.fasta\n"
	"                                             -other dir_1/               -> $PWD/dir_1/other.fasta\n"
	"                                             -other apfx                 -> $PWD/apfx.fasta\n"
	"                                             -other  (no argument)       -> aligned_out/other.fasta\n"
	"                                                                            i.e. the same output directory as used\n"
	"                                                                            for aligned output\n",
help_fastx = 
	"Output aligned reads into FASTA/FASTQ file",
help_workdir = 
	"Directory for storing Reference index,      USRDIR/sortmerna/run/\n"
	"                                            Key-value database, and the output.\n"
	"                                            Default structure:\n"
	"                                              WORKDIR/\n"
	"                                                 idx/\n"
	"                                                 kvdb/\n"
	"                                                 out/\n",
help_kvdb =
	"Directory for storing Key-value database    WORKDIR/kvdb\n"
	"                                            KVDB is used for storing alignement results.\n",
help_idx =
	"Directory for storing Reference index.      WORKDIR/idx\n"
	"                                            \n",
help_sam = 
	"Output SAM alignment for aligned reads.\n",
help_SQ = 
	"Add SQ tags to the SAM file\n",
help_blast = 
	"output alignments in various Blast-like formats\n"
	"                                            '0'                    - pairwise\n"
	"                                            '1'                    - tabular(Blast - m 8 format)\n"
	"                                            '1 cigar'              - tabular + column for CIGAR\n"
	"                                            '1 cigar qcov'         - tabular + columns for CIGAR\n"
	"                                                                     and query coverage\n"
	"                                            '1 cigar qcov qstrand' - tabular + columns for CIGAR,\n"
	"                                                                     query coverage and strand\n",
help_dbg_put_db = 
	"",
help_log = 
	"Output overall statistics.                              True\n"
	"                                            TODO: remove\n",
help_num_alignments = 
	"Positive integer (INT >=0).\n"
	"                                            Report first INT alignments per read reaching E-value threshold\n"
	"                                            If INT = 0, all alignments will be output\n"
	//"                                            Mutually exclusive with option '"+ OPT_BEST +"'\n"
	//"                                            Mutually exclusive with option '" + OPT_OTU_MAP + "'\n"
	"                                            This option allows to lower the CPU time and memory use.",
help_no_best = 
	"Disable best alignments search                          1\n"
	"                                            by searching --min_lis INT candidate alignments\n"
	"                                            If INT == 0: search All candidate alignments\n"
	"                                            If INT > 0: search INT best alignments.\n"
	"                                            The larger is the INT, the longer is the search time.\n"
	"                                            Explanation:\n"
	"                                            A read can potentially be aligned (reaching E-value threshold)\n"
	"                                            to multiple reference sequences.\n"
	"                                            The 'best' alignment is the highest scoring alignment out of All\n"
	"                                            alignments of a Read.\n"
	"                                            To find the Best alignment - an exhaustive search over All\n"
	"                                            references has to be performed.\n"
	"                                            'best 1' and 'best 0' (all the bests) are Equally intensive processes\n"
	"                                            requiring the exhaustive search. Only the size of reports will differ.\n",
help_min_lis = 
	"Search all alignments having the first INT longest LIS  2\n"
	"                                            LIS stands for Longest Increasing Subsequence,\n"
	"                                            it is computed using seeds' positions to expand hits into\n"
	"                                            longer matches prior to Smith - Waterman alignment.\n",
	//"                                            Requires option '"+ OPT_BEST +"'.\n"
	//"                                            Mutually exclusive with option '"+ OPT_NUM_ALIGNMENTS +"'\n",
help_print_all_reads = 
	"Output null alignment strings for non-aligned reads     False\n"
	"                                            to SAM and/or BLAST tabular files\n",
help_paired =
	"Input is a single file with interleaved paired reads    False\n"
	"                                            If such a file is used, and\n"
	"                                            neither '" + OPT_PAIRED_IN + "' nor '" + OPT_PAIRED_OUT + "' are specified,\n"
	"                                            use this option together with '" + OPT_OUT2 + "' to output\n"
	"                                            FWD and REV reads into separate files\n",
help_paired_in = 
	"If one of the paired-end reads is Aligned,              False\n"
	"                                            put both reads into Aligned FASTA/Q file\n"
	"                                            Must be used with '" + OPT_FASTX + "'.\n"
	"                                            Mutually exclusive with '" + OPT_PAIRED_OUT + "'.\n",
help_paired_out = 
	"If one of the paired-end reads is Non-aligned,          False\n"
	"                                            put both reads into Non-Aligned FASTA/Q file\n"
	"                                            Must be used with '" + OPT_FASTX + "'.\n"
	"                                            Mutually exclusive with '" + OPT_PAIRED_IN + "'.\n",
help_out2 =
	"Output paired reads into separate files.                False\n"
	"                                            Must be used with '" + OPT_FASTX + "'.\n"
	"                                            Ignored without either of '" + OPT_PAIRED_IN + "' |\n"
	"                                            '" + OPT_PAIRED_OUT + "' | '" + OPT_PAIRED + "' | two '" + OPT_READS + "'\n",
help_match = 
	"SW score (positive integer) for a match.                2\n",
help_mismatch = 
	"SW penalty (negative integer) for a mismatch.          -3\n",
help_gap_open = 
	"SW penalty (positive integer) for introducing a gap.    5\n",
help_gap_ext = 
	"SW penalty (positive integer) for extending a gap.      2\n",
help_N = 
	"SW penalty for ambiguous letters (N's) scored\n"
	"                                            as --mismatch\n",
help_F = 
	"Search only the forward strand.                         False\n",
help_R = 
	"Search only the reverse-complementary strand.           False\n",
help_e = 
	"E-value threshold.                                      1\n"
	"                                            Defines the 'statistical significance' of a local alignment.\n"
	"                                            Exponentially correllates with the Minimal Alignment Score.\n"
	"                                            Higher E-values (100, 1000, ...) cause More reads\n"
	"                                            to Pass the alignment threshold\n",
help_v = 
	"Produce verbose output when building the index          True\n",
help_id = 
	"%%id similarity threshold (the alignment                0.97\n"
	"                                            must still pass the E-value threshold).\n",
help_coverage = 
	"%%query coverage threshold (the alignment must          0.97\n"
	"                                            still pass the E-value threshold)\n",
help_denovo_otu = 
	"FASTA/FASTQ file for reads matching database < %%id     False\n"
	"                                            (set using --id) and < %%cov (set using --coverage)\n"
	"                                            (alignment must still pass the E-value threshold).\n",
help_otu_map = 
	"Output OTU map (input to QIIME's make_otu_table.py).    False\n",
help_passes = 
	"Three intervals at which to place the seed on the read  L,L/2,3\n"
	"                                            (L is the seed length)\n",
help_edges = 
	"Number (or percent if INT followed by %% sign) of       4\n"
	"                                            nucleotides to add to each edge of the read\n"
	"                                            prior to SW local alignment\n",
help_num_seeds = 
	"Number of seeds matched before searching                2\n"
	"                                            for candidate LIS\n",
help_pid = 
	"Add pid to output file names.                           False\n",
help_full_search = 
	"Search for all 0-error and 1-error seed                 False\n"
	"                                            matches in the index rather than stopping\n"
	"                                            after finding a 0-error match (<1%% gain in\n"
	"                                            sensitivity with up four-fold decrease in speed)\n",
help_h = 
	"Print help information\n",
help_version = 
	"Print SortMeRNA version number\n",
help_cmd = 
	"Launch an interactive session (command prompt)          False\n",
help_task = 
	"Processing Task:                                        4\n"
	"                                            0 - align. Only perform alignment\n"
	"                                            1 - post-processing (log writing)\n"
	"                                            2 - generate reports\n"
	"                                            3 - align and post-process\n"
	"                                            4 - all\n",
help_a = 
	"DEPRECATED in favour of '-threads'. Number of           numCores\n"
	"                                            processing threads to use.\n"
	"                                            Automatically redirects to '-threads'\n",
help_threads = 
	"Number of Processing threads to use                     numCores\n",
help_thpp = 
	"Number of Post-Processing Read:Process threads to use   1:1\n",
help_threp = 
	"Number of Report Read:Process threads to use            1:1\n",
help_tmpdir = 
	"Indexing: directory for writing temporary files when\n"
	"                                            building the reference index\n",
help_interval = 
	"Indexing: Positive integer: index every Nth L-mer in    1\n"
	"                                            the reference database e.g. '-interval 2'.\n",
help_m = 
	"Indexing: the amount of memory (in Mbytes) for building 3072\n"
	"                                            the index.\n",
help_L = 
	"Indexing: seed length.                                  18\n",
help_max_pos = 
	"Indexing: maximum (integer) number of positions to store  1000\n"
	"                                            for each unique L-mer. If 0 all positions are stored.\n"
;

const std::string WORKDIR_DEF_SFX = "sortmerna/run";

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
	Runopts(int argc, char** argv, bool dryrun=false);
	//~Runopts() {}

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

	// variables
public:
	// Option selection Flags
	//    alignment control
	bool is_best = true; // default if no OPT_NO_BEST was specified
	bool is_best_id_cov = false; // search for best alignments that also pass ID and COV
	bool is_min_lis = false;
	bool is_num_alignments = false; // OPT_NUM_ALIGNMENTS was specified
	bool is_full_search = false; // OPT_FULL_SEARCH was selected
	bool is_forward = false; // OPT_F was selected i.e. search only the forward strand
	bool is_reverse = false; // OPT_R was selected i.e. search only the reverse-complementary strand
	//    output control
	bool is_paired_in = false; // OPT_PAIRED_IN was selected i.e. both paired-end reads go in 'aligned' fasta/q file. Only Fasta/q and De-novo reporting.
	bool is_paired_out = false; // '--paired_out' both paired-end reads go in 'other' fasta/q file. Only Fasta/q and De-novo reporting.
	bool is_out2 = false; // 20200127 output paired reads into separate files. Issue 202
	bool is_denovo_otu = false; // output file with reads matching database < %%id (set using --id) and < %%cov (set using --coverage)
	bool is_log = true; // OPT_LOG was selected i.e. output overall statistics. TODO: remove this option, always generate.
	bool is_print_all_reads = false; // '--print_all_reads' output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files
	bool is_sam = false; // OPT_SAM was specified. output SAM alignment (for aligned reads only)
	bool is_SQ = false; // OPT_SQ add SQ tags to the SAM file
	bool is_blast = false; // OPT_BLAST was specified
	bool is_fastx = false; // OPT_FASTX was selected i.e. output FASTA/FASTQ file (for aligned and/or rejected reads)
	bool is_other = false; // OPT_OTHER was selected i.e. flags to produce 'other' file
	bool is_otu_map = false; // OPT_OTU_MAP was selected i.e. output OTU map (input to QIIME's make_otu_table.py)
	bool is_verbose; // OPT_V was selected (indexing)
	bool is_pid = false; // add pid to output file names
	bool is_cmd = false; // start interactive session
	bool is_dbg_put_kvdb = false; // if True - do Not put records into Key-value DB. Debugging Memory Consumption.

	// Option derived Flags
	bool is_as_percent = false; // derived from OPT_EDGES

	// Other flags
	bool exit_early = false; // TODO: has no action? Flag to exit processing when either the reads or the reference file is empty or not FASTA/FASTQ
	bool is_index_built = false; // flags the index is built and ready for use. TODO: this is no Option flag. Move to a more appropriate place.
	bool is_gz = false; // flags reads file is compressed and can be read. TODO: no Option related flag. Move to a proper place.
	bool is_paired = false; // flags the reads are paired

	std::filesystem::path workdir; // Directory for index, KVDB, Output
	std::filesystem::path idxdir;
	std::filesystem::path kvdbdir;
	std::filesystem::path outdir;
	std::filesystem::path aligned_pfx; // aligned reads output file prefix [dir/][pfx]
	std::filesystem::path other_pfx; // non-aligned reads output file prefix [dir/][pfx]
	std::string cmdline;

	int num_read_thread = 1; // number of threads reading the Reads file.
	int num_write_thread = 1; // number of threads writing to Key-value database
	int num_proc_thread = 0; // '-a' number of threads to use for alignment, post-processing, reporting. Default - all available cores.
	int num_read_thread_pp = 1; // number of post-processing read threads
	int num_proc_thread_pp = 1; // number of post-processing processor threads
	int num_read_thread_rep = 1; // number of report reader threads
	int num_proc_thread_rep = 1; // number of report processor threads

	int queue_size_max = 1000; // max number of Reads in the Read and Write queues. 10 works OK.

	int32_t num_alignments = 1; // [3] help_num_alignments
	int32_t hit_seeds = 2; // OPT_NUM_SEEDS Min number of seeds on a read that have matches in DB prior calculating LIS
	int32_t min_lis = 2; // OPT_MIN_LIS search all alignments having the first N longest LIS
	int32_t edges = -1; // OPT_EDGES

	uint32_t minoccur = 0; // TODO: add to cmd options. Min number of k-mer occurrences in the DB to use for matching. See 'index.lookup_tbl[kmer_idx].count'

	int match = 2; // '--match' SW score (positive integer) for a match
	int mismatch = -3; // '--mismatch' SW penalty (negative integer) for a mismatch
	long gap_open = 5; // '--gap_open' SW penalty (positive integer) for introducing a gap
	long gap_extension = 2; // '--gap_ext' SW penalty (positive integer) for extending a gap
	int score_N = 0; // '-N' SW penalty for ambiguous letters (N's)

	double evalue = -1.0; // '-e' E-value threshold
	double min_id = -1.0; // OTU-picking option: Identity threshold (%ID)
	double min_cov = -1.0; // [4] OTU-picking option: minimum Coverage (%COV)

	// indexing options
	double max_file_size = 3072; // max size of an index file (or a part of the file). When exceeded, the index is split into parts.
	uint32_t seed_win_len = 18; // OPT_L seed kmer length
	uint32_t interval = 1; // size of k-mer window shift. Default 1 is the min possible to generate max number of k-mers.
	uint32_t max_pos = 10000;
	// ~ END indexing options

	std::vector<std::string> blastops; // [1]
	std::vector<std::string> readfiles; // '--reads'
	std::vector<std::pair<std::string, std::string>> indexfiles; // '-ref' pairs 'Ref_file:Idx_file_pfx'
	std::vector<std::vector<uint32_t>> skiplengths; // [2] OPT_PASSES K-mer window shift sizes. Refstats::load

	const std::string dbkey = "run_options";
	const std::string IDX_DIR  = "idx";
	const std::string KVDB_DIR = "kvdb";
	const std::string OUT_DIR  = "out";

	enum ALIGN_REPORT { align, postproc, report, alipost, all };
	ALIGN_REPORT alirep = ALIGN_REPORT::all;
	BlastFormat blastFormat = BlastFormat::TABULAR;

	// methods
private:
	void process(int argc, char**argv, bool dryrun);
	void validate();
	void validate_idxdir(); // called from validate
	void validate_kvdbdir(); // called from validate
	void validate_aligned_pfx();
	void validate_other_pfx();

	void opt_sort();
	void opt_reads(const std::string& val);
	void opt_reads_gz(char **argv, int& narg);
	void opt_ref(const std::string& val);
	void opt_aligned(const std::string &val);
	void opt_other(const std::string &val);
	void opt_log(const std::string &val);
	void opt_denovo_otu(const std::string &val);
	void opt_otu_map(const std::string &val);
	void opt_print_all_reads(const std::string &val);
	void opt_pid(const std::string &val);
	void opt_paired(const std::string& val);
	void opt_paired_in(const std::string &val);
	void opt_paired_out(const std::string &val);
	void opt_out2(const std::string& val);
	void opt_match(const std::string &val);
	void opt_mismatch(const std::string &val);
	void opt_gap_open(const std::string &val);
	void opt_gap_ext(const std::string &val);
	void opt_num_seeds(const std::string &val);
	void opt_fastx(const std::string &val);
	void opt_sam(const std::string& val);
	void opt_blast(const std::string& val);
	void opt_min_lis(const std::string& val);
	void opt_no_best(const std::string& val);
	void opt_num_alignments(const std::string& val);
	void opt_edges(const std::string& val);
	void opt_full_search(const std::string& val);
	void opt_SQ(const std::string& val);
	void opt_passes(const std::string& val);
	void opt_id(const std::string& val);
	void opt_coverage(const std::string& val);
	void opt_version(const std::string& val);
	void opt_task(const std::string& val);
	void opt_cmd(const std::string& val);
	void opt_threads(const std::string& val);
	void opt_thpp(const std::string& val); // post-proc threads --thpp 1:1
	void opt_threp(const std::string& val); // report threads --threp 1:1 
	void opt_a(const std::string& val);
	void opt_e(const std::string& val); // opt_e_Evalue
	void opt_F(const std::string& val); // opt_F_ForwardOnly
	void opt_R(const std::string& val); // opt_R_ReverseOnly
	void opt_h(const std::string& val);
	void opt_v(const std::string& val); // opt_v_Verbose
	void opt_N(const std::string& val); // opt_N_MatchAmbiguous
	void opt_workdir(const std::string& path);
	void opt_kvdb(const std::string& path);
	void opt_idx(const std::string& path);

	// ref tmpdir interval m L max_pos v h  // indexing options
	void opt_tmpdir(const std::string &val);
	void opt_interval(const std::string &val);
	void opt_m(const std::string &val);
	void opt_L(const std::string &val);
	void opt_max_pos(const std::string &val);

	void opt_default(const std::string& opt);
	void opt_dbg_put_db(const std::string& opt);
	void opt_unknown(char** argv, int& narg, char* opt);

	std::string to_string();
	std::string to_bin_string();
	void store_to_db(KeyValueDatabase& kvdb);

	// variables
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
	bool have_reads = false; // flags reads file is plain text and can be read

	// container for options passed to the program
	std::multimap<std::string, std::string> mopt;

	// OPTIONS Map - specifies all possible options
	const std::array<opt_6_tuple, 48> options = {
		std::make_tuple(OPT_REF,            "PATH",        COMMON,      true,  help_ref, &Runopts::opt_ref),
		std::make_tuple(OPT_READS,          "PATH",        COMMON,      true,  help_reads, &Runopts::opt_reads),
		std::make_tuple(OPT_WORKDIR,        "PATH",        COMMON,      false, help_workdir, &Runopts::opt_workdir),
		std::make_tuple(OPT_KVDB,           "PATH",        COMMON,      false, help_kvdb, &Runopts::opt_kvdb),
		std::make_tuple(OPT_IDX,            "PATH",        COMMON,      false, help_idx, &Runopts::opt_idx),
		std::make_tuple(OPT_FASTX,          "BOOL",        COMMON,      false, help_fastx, &Runopts::opt_fastx),
		std::make_tuple(OPT_SAM,            "BOOL",        COMMON,      false, help_sam, &Runopts::opt_sam),
		std::make_tuple(OPT_SQ,             "BOOL",        COMMON,      false, help_SQ, &Runopts::opt_SQ),
		std::make_tuple(OPT_BLAST,          "STRING",      COMMON,      false, help_blast, &Runopts::opt_blast),
		std::make_tuple(OPT_ALIGNED,        "STRING/BOOL", COMMON,      false, help_aligned, &Runopts::opt_aligned),
		std::make_tuple(OPT_OTHER,          "STRING/BOOL", COMMON,      false, help_other, &Runopts::opt_other),
		std::make_tuple(OPT_NUM_ALIGNMENTS, "INT",         COMMON,      false, help_num_alignments, &Runopts::opt_num_alignments),
		std::make_tuple(OPT_NO_BEST,        "BOOL",        COMMON,      false, help_no_best, &Runopts::opt_no_best),
		std::make_tuple(OPT_MIN_LIS,        "INT",         COMMON,      false, help_min_lis, &Runopts::opt_min_lis),
		std::make_tuple(OPT_PRINT_ALL_READS,"BOOL",        COMMON,      false, help_print_all_reads, &Runopts::opt_print_all_reads),
		std::make_tuple(OPT_PAIRED,         "BOOL",        COMMON,      false, help_paired, &Runopts::opt_paired),
		std::make_tuple(OPT_PAIRED_IN,      "BOOL",        COMMON,      false, help_paired_in, &Runopts::opt_paired_in),
		std::make_tuple(OPT_PAIRED_OUT,     "BOOL",        COMMON,      false, help_paired_out, &Runopts::opt_paired_out),
		std::make_tuple(OPT_OUT2,           "BOOL",        COMMON,      false, help_out2, &Runopts::opt_out2),
		std::make_tuple(OPT_MATCH,          "INT",         COMMON,      false, help_match, &Runopts::opt_match),
		std::make_tuple(OPT_MISMATCH,       "INT",         COMMON,      false, help_mismatch, &Runopts::opt_mismatch),
		std::make_tuple(OPT_GAP_OPEN,       "INT",         COMMON,      false, help_gap_open, &Runopts::opt_gap_open),
		std::make_tuple(OPT_GAP_EXT,        "INT",         COMMON,      false, help_gap_ext, &Runopts::opt_gap_ext),
		std::make_tuple(OPT_E,              "DOUBLE",      COMMON,      false, help_e, &Runopts::opt_e),
		std::make_tuple(OPT_F,              "BOOL",        COMMON,      false, help_F, &Runopts::opt_F),
		std::make_tuple(OPT_N,              "BOOL",        COMMON,      false, help_N, &Runopts::opt_N),
		std::make_tuple(OPT_R,              "BOOL",        COMMON,      false, help_R, &Runopts::opt_R),
		std::make_tuple(OPT_ID,             "INT",         OTU_PICKING, false, help_id, &Runopts::opt_id),
		std::make_tuple(OPT_COVERAGE,       "INT",         OTU_PICKING, false, help_coverage, &Runopts::opt_coverage),
		std::make_tuple(OPT_DENOVO_OTU,     "BOOL",        OTU_PICKING, false, help_denovo_otu, &Runopts::opt_denovo_otu),
		std::make_tuple(OPT_OTU_MAP,        "BOOL",        OTU_PICKING, false, help_otu_map, &Runopts::opt_otu_map),
		std::make_tuple(OPT_PASSES,         "INT,INT,INT", ADVANCED,    false, help_passes, &Runopts::opt_passes),
		std::make_tuple(OPT_EDGES,          "INT",         ADVANCED,    false, help_edges, &Runopts::opt_edges),
		std::make_tuple(OPT_NUM_SEEDS,      "BOOL",        ADVANCED,    false, help_num_seeds, &Runopts::opt_num_seeds),
		std::make_tuple(OPT_FULL_SEARCH,    "INT",         ADVANCED,    false, help_full_search, &Runopts::opt_full_search),
		std::make_tuple(OPT_PID,            "BOOL",        ADVANCED,    false, help_pid, &Runopts::opt_pid),
		std::make_tuple(OPT_A,              "INT",         ADVANCED,    false, help_a, &Runopts::opt_a),
		std::make_tuple(OPT_THREADS,        "INT",         ADVANCED,    false, help_threads, &Runopts::opt_threads),
		std::make_tuple(OPT_L,              "DOUBLE",      INDEXING,    false, help_L, &Runopts::opt_L),
		std::make_tuple(OPT_M,              "DOUBLE",      INDEXING,    false, help_m, &Runopts::opt_m),
		std::make_tuple(OPT_V,              "BOOL",        INDEXING,    false, help_v, &Runopts::opt_v),
		std::make_tuple(OPT_INTERVAL,       "INT",         INDEXING,    false, help_interval, &Runopts::opt_interval),
		std::make_tuple(OPT_MAX_POS,        "INT",         INDEXING,    false, help_max_pos, &Runopts::opt_max_pos),
		std::make_tuple(OPT_H,              "BOOL",        HELP,        false, help_h, &Runopts::opt_h),
		std::make_tuple(OPT_VERSION,        "BOOL",        HELP,        false, help_version, &Runopts::opt_version),
		std::make_tuple(OPT_DBG_PUT_DB,     "BOOL",        DEVELOPER,   false, help_dbg_put_db, &Runopts::opt_dbg_put_db),
		std::make_tuple(OPT_CMD,            "BOOL",        DEVELOPER,   false, help_cmd, &Runopts::opt_cmd),
		std::make_tuple(OPT_TASK,           "INT",         DEVELOPER,   false, help_task, &Runopts::opt_task)
		//std::make_tuple(OPT_THPP,           "INT:INT",     DEVELOPER,   false, help_thpp, &Runopts::opt_thpp),
		//std::make_tuple(OPT_THREP,          "INT:INT",     DEVELOPER,   false, help_threp, &Runopts::opt_threp)
	};
	// ~map options
}; // ~struct Runopts
// ~options.cpp
