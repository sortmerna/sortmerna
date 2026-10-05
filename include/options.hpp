/*
@copyright 2016-2026 Clarity Genomics Inc
@copyright 2012-2016 Bonsai Bioinformatics Research Group
@copyright 2014-2016 Knight Lab, Department of Pediatrics, UCSD, La Jolla

@parblock
SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA

This is a free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SortMeRNA is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
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

/*
 * file: options.hpp
 * created: Aug 19, 2017 Sat
 *
 * skiplength
 *    skip lengths for pass 1, pass 2 and pass 3 in first step of sortmerna
 *    pipeline for each reference database searched
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <array>
#include <filesystem>
#include <cstdint>  // uint64_t

#include "common.hpp"

namespace sortmerna {

/*! @brief Maximum length of input reads
    (not limited to this length algorithmically)
*/
const uint64_t MAX_READ_LEN = 30000;

/* Canonical option names.
 *
 * Canonical spelling is dash-separated. Underscore spellings of the same name
 * are accepted on the command line and resolve to the canonical form (see
 * canonicalize_opt_name/normalize_argv in options.cpp), but are not advertised
 * in the help. Single-character names are short options i.e. '-N', not '--N'.
 */
const std::string \
OPT_REF = "ref",
OPT_READS = "reads",
OPT_ALIGNED = "aligned",
OPT_OTHER = "other",
OPT_WORKDIR = "workdir",
OPT_KVDB = "kvdb",
OPT_IDXDIR = "idx-dir",
OPT_READB = "readb",
OPT_FASTX = "fastx",
OPT_SAM = "sam",
OPT_SQ = "SQ",
OPT_BLAST = "blast",
OPT_NUM_ALIGNMENTS = "num-alignments",
OPT_NO_BEST = "no-best",
OPT_MIN_LIS = "min-lis",
OPT_PRINT_ALL_READS = "print-all-reads",
OPT_PAIRED = "paired",
OPT_PAIRED_IN = "paired-in",
OPT_PAIRED_OUT = "paired-out",
OPT_OUT2 = "out2",
OPT_SOUT = "sout",
OPT_MATCH = "match",
OPT_MISMATCH = "mismatch",
OPT_GAP_OPEN = "gap-open",
OPT_GAP_EXT = "gap-ext",
OPT_A = "a",
OPT_E = "e",
OPT_F = "F",
OPT_H = "h",
OPT_HELP = "help",
OPT_L = "L",
OPT_M = "m",
OPT_N = "N",
OPT_R = "R",
OPT_V = "v",
OPT_ID = "id",
OPT_COVERAGE = "coverage",
OPT_DENOVO_OTU = "de-novo-otu",
OPT_OTU_MAP = "otu-map",
OPT_PASSES = "passes",
OPT_EDGES = "edges",
OPT_NUM_SEEDS = "num-seeds",
OPT_FULL_SEARCH = "full-search",
OPT_PID = "pid",
OPT_VERSION = "version",
OPT_CMD = "cmd",
OPT_TASK = "task",
OPT_THREADS = "threads",
OPT_DBG_PUT_DB = "dbg-put-db",
OPT_INTERVAL = "interval",
OPT_MAX_POS = "max-pos",
OPT_READFEED = "readfeed",
OPT_ZIP_OUT = "zip-out",
OPT_INDEX = "index",
OPT_ALIGN = "align",  // TODO: on hold
OPT_FILTER = "filter",  // TODO: on hold
OPT_DBG_LEVEL = "dbg-level",
OPT_MAX_READ_LEN = "max-read-len",
OPT_SCORE_SPLIT = "score-split",
OPT_FLUSH_DELAY = "flush-delay";

/* Help text.
 *
 * Two strings per option:
 *   help_<opt>        one line, no padding, no trailing newline. Used by '-h'.
 *   help_<opt>_long   optional extended body. Used by '--help' and by the
 *                     per-option help ('sortmerna <option> -h').
 *
 * Default values are NOT written into these strings - they live in
 * OptDesc::def so the renderer can place them in their own column.
 */
const std::string \
help_header =
"  Usage:   sortmerna -ref FILE [-ref FILE] -reads FWD_READS [-reads REV_READS] [OPTIONS]\n"
"\n"
"  Option names may be given with one or two leading dashes ('-ref' == '--ref'),\n"
"  and dashes and underscores are interchangeable ('--paired-out' == '--paired_out').\n"
"  A value may follow the name or be attached with '=' ('--threads 4' == '--threads=4').\n"
"\n"
"  '-h' lists the options, '--help' prints this full reference, and\n"
"  'sortmerna <option> -h' explains a single option e.g. 'sortmerna index -h'.\n",

help_ref =
	"Reference file (FASTA) absolute or relative path.",
help_ref_long =
	"Use multiple times, once per a reference file.",

help_reads =
	"Raw reads file (FASTA/FASTQ/FASTA.GZ/FASTQ.GZ).",
help_reads_long =
	"Use twice for files with paired reads.\n"
	"The file extensions are Not important. The program automatically\n"
	"recognizes the file format as flat/compressed, fasta/fastq.",

help_aligned =
	"Aligned reads file prefix [dir/][pfx]",
help_aligned_long =
	"Directory and file prefix for aligned output i.e. each\n"
	"output file goes into the specified directory with the given prefix.\n"
	"The appropriate extension: (fasta|fastq|blast|sam|etc) is automatically added.\n"
	"Both 'dir' and 'pfx' are optional.\n"
	"The 'dir' can be a relative or an absolute path.\n"
	"If 'dir' is not specified, the output is created in the WORKDIR/out/\n"
	"If 'pfx' is not specified, the prefix 'aligned' is used\n"
	"Examples:\n"
	"  '-aligned $MYDIR/dir_1/dir_2/1' -> $MYDIR/dir_1/dir_2/1.fasta\n"
	"  '-aligned dir_1/apfx'           -> $PWD/dir_1/apfx.fasta\n"
	"  '-aligned dir_1/'               -> $PWD/aligned.fasta\n"
	"  '-aligned apfx'                 -> $PWD/apfx.fasta\n"
	"  '-aligned  (no argument)'       -> WORKDIR/out/aligned.fasta",

help_other =
	"Non-aligned reads file prefix [dir/][pfx]",
help_other_long =
	"Directory and file prefix for non-aligned output i.e. each\n"
	"output file goes into the specified directory with the given prefix.\n"
	"The appropriate extension: (fasta|fastq|blast|sam|etc) is automatically added.\n"
	"Must be used with '" + OPT_FASTX + "'.\n"
	"Both 'dir' and 'pfx' are optional.\n"
	"The 'dir' can be a relative or an absolute path.\n"
	"If 'dir' is not specified, the output is created in the WORKDIR/out/\n"
	"If 'pfx' is not specified, the prefix 'other' is used\n"
	"Examples:\n"
	"  '-other $MYDIR/dir_1/dir_2/1' -> $MYDIR/dir_1/dir_2/1.fasta\n"
	"  '-other dir_1/apfx'           -> $PWD/dir_1/apfx.fasta\n"
	"  '-other dir_1/'               -> $PWD/dir_1/other.fasta\n"
	"  '-other apfx'                 -> $PWD/apfx.fasta\n"
	"  '-other  (no argument)'       -> aligned_out/other.fasta\n"
	"                                   i.e. the same output directory\n"
	"                                   as used for aligned output",

help_fastx =
	"Output aligned reads into FASTA/FASTQ file",

help_workdir =
	"Workspace directory",
help_workdir_long =
	"Default structure: WORKDIR/\n"
	"                       idx/   (References index)\n"
	"                       kvdb/  (Key-value storage for alignments)\n"
	"                       out/   (processing output)\n"
	"                       readb/ (pre-processed reads/index)",

help_kvdb =
	"Directory for Key-value database",
help_kvdb_long =
	"KVDB is used for storing the alignment results.",

help_idxdir =
	"Directory for storing Reference index.",

help_readb =
	"Storage for pre-processed reads",
help_readb_long =
	"Directory storing the split reads.\n"
	"Only used with '-" + OPT_READFEED + "' = "
	    + std::to_string(static_cast<int>(FEED_TYPE::SPLIT_READS)) +
	    " i.e. Split Reads. Otherwise ignored.",

help_readfeed =
	"Method of accessing the reads by the reads processors",
help_readfeed_long =
	"0 - Indexed. Random access to compressed or flat reads files using\n"
	"    pre-calculated indices\n"
	"1 - DEPRECATED: Split reads. Reads files are split into parts equal\n"
	"    the number of processing threads",

help_sam =
	"Output SAM alignment for aligned reads.",
help_SQ =
	"Add SQ tags to the SAM file",

help_blast =
	"Output alignments in various Blast-like formats",
help_blast_long =
	"Sample values: '0'                    - pairwise\n"
	"               '1'                    - tabular (Blast -m 8 format)\n"
	"               '1 cigar'              - tabular + column for CIGAR\n"
	"               '1 cigar qcov'         - tabular + columns for CIGAR and query coverage\n"
	"               '1 cigar qcov qstrand' - tabular + columns for CIGAR, query coverage,\n"
	"                                        and strand",

help_dbg_put_db =
	"Do Not put records into the Key-value DB (debugging memory consumption)",

help_num_alignments =
	"Number of alignments to output per read (INT >= 0).",
help_num_alignments_long =
	"If used with '-" + OPT_NO_BEST + "' reports first INT alignments per read reaching\n"
	"E-value threshold, which allows to lower the CPU time and memory use.\n"
	"Otherwise outputs INT best alignments.\n"
	"If INT = 0, all alignments are output.",

help_no_best =
	"Disable best alignments search",
help_no_best_long =
	"The 'best' alignment is the highest scoring alignment out of All alignments of a read,\n"
	"and the read can potentially be aligned (reaching E-value threshold) to multiple reference\n"
	"sequences.\n"
	"By default the program searches for best alignments i.e. performs an exhaustive search\n"
	"over all references. Using '-" + OPT_NO_BEST + "' will make the program to search just\n"
	"the first N alignments, where N is set using '-" + OPT_NUM_ALIGNMENTS + "' i.e. 1 by default.",

help_min_lis =
	"Search only alignments that have the LIS of at least N seeds long",
help_min_lis_long =
	"LIS stands for Longest Increasing Subsequence. It is computed using seeds, which\n"
	"are k-mers common to the read and the reference sequence. Sorted sequences of such seeds\n"
	"are used to filter the candidate references prior performing the Smith-Waterman alignment.",

help_print_all_reads =
	"Output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files",

help_paired =
	"Flags paired reads",
help_paired_long =
	"If a single reads file is provided, use this option to indicate\n"
	"the file contains interleaved paired reads when neither\n"
	"'" + OPT_PAIRED_IN + "' | '" + OPT_PAIRED_OUT + "' | '" + OPT_OUT2 + "' | '" + OPT_SOUT + "' are specified.",

help_paired_in =
	"Flags the paired-end reads as Aligned, when either of them is Aligned.",
help_paired_in_long =
	"With this option both reads are output into Aligned FASTA/Q file.\n"
	"Must be used with '" + OPT_FASTX + "'.\n"
	"Mutually exclusive with '" + OPT_PAIRED_OUT + "'.",

help_paired_out =
	"Flags the paired-end reads as Non-aligned, when either of them is non-aligned.",
help_paired_out_long =
	"With this option both reads are output into Non-Aligned FASTA/Q file.\n"
	"Must be used with '" + OPT_FASTX + "'.\n"
	"Mutually exclusive with '" + OPT_PAIRED_IN + "'.",

help_out2 =
	"Output paired reads into separate files.",
help_out2_long =
	"Must be used with '" + OPT_FASTX + "'.\n"
	"If a single reads file is provided, this option implies interleaved paired reads.\n"
	"When used with '" + OPT_SOUT + "', four (4) output files for aligned reads will be generated:\n"
	"'aligned-paired-fwd, aligned-paired-rev, aligned-singleton-fwd, aligned-singleton-rev'.\n"
	"If '" + OPT_OTHER + "' option is also used, eight (8) output files will be generated.",

help_sout =
	"Separate paired and singleton aligned reads.",
help_sout_long =
	"To be used with '" + OPT_FASTX + "'.\n"
	"If a single reads file is provided, this option implies interleaved paired reads.\n"
	"Cannot be used with '" + OPT_PAIRED_IN + "' | '" + OPT_PAIRED_OUT + "'.",

help_match =
	"SW score (positive integer) for a match.",
help_mismatch =
	"SW penalty (negative integer) for a mismatch.",
help_gap_open =
	"SW penalty (positive integer) for introducing a gap.",
help_gap_ext =
	"SW penalty (positive integer) for extending a gap.",
help_N =
	"SW penalty for ambiguous letters (N's) scored as '--" + OPT_MISMATCH + "'.",
help_F =
	"Search only the forward strand.",
help_R =
	"Search only the reverse-complementary strand.",

help_e =
	"E-value threshold.",
help_e_long =
	"Defines the 'statistical significance' of a local alignment.\n"
	"Exponentially correlates with the Minimal Alignment score.\n"
	"Higher E-values (100, 1000, ...) cause More reads to Pass the alignment threshold.",

help_v =
	"Produce verbose output when building the index",

help_id =
	"%id similarity threshold (the alignment must still pass the E-value threshold).",
help_coverage =
	"%query coverage threshold (the alignment must still pass the E-value threshold)",

help_denovo_otu =
	"Output FASTA file with 'de novo' reads",
help_denovo_otu_long =
	"Read is 'de novo' if its alignment score passes E-value threshold, but both the identity\n"
	"'-" + OPT_ID + "', and the '-" + OPT_COVERAGE + "' are below their corresponding thresholds\n"
	"i.e. ID < %id and COV < %cov.",

help_otu_map =
	"Output OTU map (input to QIIME's make_otu_table.py).",
help_otu_map_long =
	"Cannot be used with '" + OPT_NO_BEST + "' because the grouping is done around\n"
	"the best alignment.",

help_passes =
	"Three intervals at which to place the seed on the read (L is the seed length)",

help_edges =
	"Number (or percent if INT followed by % sign) of nucleotides to add to each edge",
help_edges_long =
	"The nucleotides are added to each edge of the read prior to SW local alignment.",

help_num_seeds =
	"Number of seeds matched before searching for candidate LIS",

help_pid =
	"Add pid to output file names.",

help_full_search =
	"Search for all 0-error and 1-error seed matches in the index",
help_full_search_long =
	"Rather than stopping after finding a 0-error match (<1% gain in\n"
	"sensitivity with up four-fold decrease in speed).",

help_h =
	"Print the list of options. Use 'sortmerna <option> -h' for help on one option",
help_help =
	"Print the full help reference",
help_version =
	"Print SortMeRNA version number",
help_cmd =
	"Launch an interactive session (command prompt)",

help_task =
	"Processing Task",
help_task_long =
	"Possible values: 0 - align. Only perform alignment\n"
	"                 1 - post-processing (log writing)\n"
	"                 2 - generate reports\n"
	"                 3 - align and post-process\n"
	"                 4 - all\n"
	"                 5 - index only",

help_a =
	"DEPRECATED in favour of '--" + OPT_THREADS + "'. Number of processing threads to use.",
help_a_long =
	"Automatically redirects to '--" + OPT_THREADS + "'.",

help_threads =
	"Number of Processing threads to use",

help_interval =
	"Indexing: index every Nth L-mer in the reference database e.g. '-interval 2'.",
help_m =
	"Indexing: the amount of memory (in Mbytes) for building the index.",
help_L =
	"Indexing: seed length.",

help_max_pos =
	"Indexing: maximum number of positions to store for each unique L-mer.",
help_max_pos_long =
	"If 0 - all positions are stored.",

help_zip_out =
	"Controls the output compression",
help_zip_out_long =
	"By default the report files are produced in the same format as the input i.e.\n"
	"if the reads files are compressed (gz), the output is also compressed.\n"
	"The default behaviour can be overridden by using '-" + OPT_ZIP_OUT + "'.\n"
	"The possible values: '1/true/t/yes/y'\n"
	"                     '0/false/f/no/n'\n"
	"                     '-1' (the same format as input - default)\n"
	"The values are Not case sensitive i.e. 'Yes, YES, yEs, Y, y' are all OK.\n"
	"Examples:\n"
	"  '-" + OPT_READS + " freads.gz -" + OPT_ZIP_OUT + " n' : generate flat output when the input is compressed\n"
	"  '-" + OPT_READS + " freads.flat -" + OPT_ZIP_OUT + "' : compress the output when the input files are flat",

help_index =
	"Build reference database index",
help_index_long =
	"By default when this option is not used, the program checks the reference index and\n"
	"builds it if not already existing.\n"
	"This can be changed by using '-" + OPT_INDEX + "' as follows:\n"
	"  '-" + OPT_INDEX + " 0' - skip indexing. If the index does not exist, the program will\n"
	"               terminate and warn to build the index prior performing the alignment\n"
	"  '-" + OPT_INDEX + " 1' - the default behaviour, the same as when not using this option at all",

help_dbg_level =
	"Debug level",
help_dbg_level_long =
	"Controls verbosity of the execution trace. Default value of 0 corresponds to\n"
	"the least verbose output. The highest value currently is 2.",

help_max_read_len =
	"Maximum allowed read length",

help_score_split =
	"Calculate minimal SW score per split rather than all reads.",
help_score_split_long =
	"This has an effect similar to increasing e-value i.e. lowers the filtering\n"
	"threshold to less sensitive (see issue 453).",

help_flush_delay =
	"Mid-alignment restart checkpoint interval (seconds).",
help_flush_delay_long =
	"Interval in seconds to flush/persist the alignment related data to the kvdb,\n"
	"so that the alignment could resume from the point of interrupt rather than\n"
	"repeating the process from the start."
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

	/* Does the option take a value on the command line?
	 *   NONE     - a flag. The handler ignores the value.
	 *   REQUIRED - the handler needs a value.
	 *   OPTIONAL - the handler works either way e.g. '-aligned' vs '-aligned dir/pfx'.
	 */
	enum class OptArg { NONE, REQUIRED, OPTIONAL };

	typedef void (Runopts::*OptsMemFunc)(const std::string&); // pointer to member function

	struct OptDesc
	{
		std::string  name;        // canonical name, dash-separated
		std::string  type;        // value type shown in the help e.g. "PATH", "INT"
		OPT_CATEGORY cat;         // help category
		bool         required;    // the program cannot run without it
		OptArg       arg;         // whether a value is expected
		std::string  def;         // default value, for display only
		std::string  help;        // one line summary
		std::string  help_long;   // extended description, may be empty
		OptsMemFunc  fn;          // option processing function
	};

	const std::map<OPT_CATEGORY, std::string> opt_category_name_map {
		{COMMON, "COMMON"},
		{OTU_PICKING, "OTU_PICKING"},
		{ADVANCED, "ADVANCED"},
		{DEVELOPER, "DEVELOPER"},
		{HELP, "HELP"},
		{INDEXING, "INDEXING"}
	};

	/* print the compact option list ('-h') */
	void print_help();
	/* print the full help reference ('--help') */
	void print_help_full();
	/* print the help for a single option ('sortmerna <option> -h') */
	void print_option_help(const OptDesc& opt);

	// variables
public:
	// Option selection Flags
	//    alignment control
	bool is_best = true; // default if no OPT_NO_BEST was specified
	bool is_best_id_cov = false; // TODO: search for best alignments that also pass ID and COV. Not yet implemented 20200703
	bool is_min_lis = false;
	bool is_num_alignments = false; // OPT_NUM_ALIGNMENTS was specified
	bool is_full_search = false; // OPT_FULL_SEARCH was selected
	bool is_forward = false; // OPT_F was selected i.e. search only the forward strand
	bool is_reverse = false; // OPT_R was selected i.e. search only the reverse-complementary strand
	//    output control
	bool is_paired_in = false; // OPT_PAIRED_IN was selected i.e. both paired-end reads go in 'aligned' fasta/q file. Only Fasta/q and De-novo reporting.
	bool is_paired_out = false; // '--paired-out' both paired-end reads go in 'other' fasta/q file. Only Fasta/q and De-novo reporting.
	bool is_out2 = false; // 20200127 output paired reads into separate files. Issue 202
	bool is_sout = false; // 20210105 separate singletons and paired
	bool is_otu_map = false; // OPT_OTU_MAP was selected i.e. output OTU map (input to QIIME's make_otu_table.py)
	bool is_denovo = false; // output file with reads matching database < %id (set using --id) and < %cov (set using --coverage)
	bool is_log = true; // output overall statistics. Always on.
	bool is_print_all_reads = false; // '--print-all-reads' output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files
	bool is_sam = false; // OPT_SAM was specified. output SAM alignment (for aligned reads only)
	bool is_SQ = false; // OPT_SQ add SQ tags to the SAM file
	bool is_blast = false; // OPT_BLAST was specified
	bool is_fastx = false; // OPT_FASTX was selected i.e. output FASTA/FASTQ file (for aligned and/or rejected reads)
	bool is_other = false; // OPT_OTHER was selected i.e. flags to produce 'other' file
	bool is_verbose; // OPT_V was selected (indexing)
	bool is_pid = false; // add pid to output file names
	bool is_cmd = false; // start interactive session
	bool is_dbg_put_kvdb = false; // if True - do Not put records into Key-value DB. Debugging Memory Consumption.
	int  findex = 1; // 0 (don't build index) | 1 (default - build index if not present)
	bool is_align = false;
	bool is_filter = false;
    bool is_score_split = false;  // if true - calculate the SW score per split rather then for all reads
    /* Set by the C API (smr_api) only; there is no CLI option. The library
     * aligns many independent batches against one loaded index and needs the
     * per-read result to be the same however the caller batches its reads:
     *  - Readstats neither restores nor stores its counts in the kvdb
     *  - Refstats sets minimal_score to 0 (the SW-score threshold depends on
     *    the batch's total read length); callers filter on e-value instead
     *  - align2 does not restore reads from the kvdb, always stores each
     *    read's result, and writes no restart progress */
    bool is_library_mode = false;

	// Option derived Flags
	bool is_as_percent = false; // derived from OPT_EDGES

	// Other flags
	bool exit_early = false; // TODO: has no action? Flag to exit processing when either the reads or the reference file is empty or not FASTA/FASTQ
	bool is_index_built = false; // flags the index is built and ready for use. TODO: this is no Option flag. Move to a more appropriate place.
	//bool is_gz = false; // flags reads file is compressed and can be read. TODO: no Option related flag. Move to a proper place.
	bool is_paired = false; // flags the reads are paired

	std::filesystem::path workdir; // Directory for index, KVDB, Output
	std::filesystem::path idxdir;
	std::filesystem::path kvdbdir;
	std::filesystem::path outdir;
	std::filesystem::path readb_dir; // Reads DB directory to use for split reads or split reads index. See option REED_FEED.
	std::filesystem::path aligned_pfx; // aligned reads output file prefix [dir/][pfx]
	std::filesystem::path other_pfx; // non-aligned reads output file prefix [dir/][pfx]
	std::string cmdline;

	unsigned num_read_thread = 1;  // number of threads reading the Reads file.
	unsigned num_write_thread = 1;  // number of threads writing to Key-value database
	unsigned num_proc_thread = 2;  // number of threads to use for alignment
	unsigned num_read_thread_pp = 1;  // number of post-processing read threads
	unsigned num_proc_thread_pp = 1;  // number of post-processing processor threads
	unsigned num_read_thread_rep = 1;  // number of report reader threads
	unsigned num_proc_thread_rep = 1;  // number of report processor threads
	unsigned dbg_level = 0;  // lowest debug level - minimal info.
	unsigned flush_delay = 60;  // seconds between mid-pass restart watermark commits. See restart.hpp.

	unsigned queue_size_max = 1000; // max number of Reads in the Read and Write queues. 10 works OK.
    uint64_t max_read_len = MAX_READ_LEN; // max allowed read len
	/*
	* 0 (false) | 1 (true) | -1 (not set)
	* read.is_zip  zip_out  out_zip
	* -----------------------------
	*      1         -1       1    zip
	*      1          0       0    flat
	*      1          1       1    zip
	*      0         -1       0    flat
	*      0          0       0    flat
	*      0          1       1    zip
	*/
	int zip_out = -1;

	uint32_t num_alignments = 1; // [3] help_num_alignments
	int32_t num_seeds = 2; // min number of seeds on a read that have matches in DB prior calculating LIS
	int32_t min_lis = 2; // search all alignments that have LIS >= min_lis
	int32_t edges = -1; // OPT_EDGES

	uint32_t minoccur = 0; // TODO: add to cmd options. Min number of k-mer occurrences in the DB to use for matching. See 'index.lookup_tbl[kmer_idx].count'

	int match = 2; // '--match' SW score (positive integer) for a match
	int mismatch = -3; // '--mismatch' SW penalty (negative integer) for a mismatch
	long gap_open = 5; // '--gap-open' SW penalty (positive integer) for introducing a gap
	long gap_extension = 2; // '--gap-ext' SW penalty (positive integer) for extending a gap
	int score_N = 0; // '-N' SW penalty for ambiguous letters (N's)
	FEED_TYPE feed_type = FEED_TYPE::INDEXED; // OPT_READFEED

	double evalue = -1.0; // '-e' E-value threshold
	double min_id = -1.0; // OTU-picking option: Identity threshold (%ID)
	double min_cov = -1.0; // [4] OTU-picking option: minimum Coverage (%COV)

	// indexing options
	double max_file_size = 3072; // max size of an index file in MB (or a part of the file). When exceeded, the index is split into parts.
	uint32_t seed_win_len = 18; // OPT_L seed lmer length
	uint32_t interval = 1; // size of k-mer window shift. Default 1 is the min possible to generate max number of k-mers.
	uint32_t max_pos = 10000;
	// ~ END indexing options

	std::vector<std::string> blastops; // [1]
	std::vector<std::string> readfiles; // '--reads'
	// list of pairs<ref_file, idx_file_pfx>
	//                 |         |_populated during indexing
	//                 |_populated during options processing
	std::vector<std::pair<std::string, std::string>> indexfiles;
	std::vector<std::vector<uint32_t>> skiplengths; // [2] OPT_PASSES K-mer window shift sizes. Refstats::load
	// The 3 skip lengths parsed from '--passes'. Stored here by opt_passes and
	// replicated into 'skiplengths' (per index) at the end of process(), because
	// options are processed in alphabetical order ('passes' before 'ref'), so
	// 'indexfiles' is still empty while opt_passes runs.
	std::vector<uint32_t> skiplengths_user;

	const std::string dbkey = "run_options";
	const std::string IDX_DIR  = "idx";
	const std::string KVDB_DIR = "kvdb";
	const std::string OUT_DIR  = "out";
	const std::string READB_DIR = "readb";

	enum TASK { align, summary, report, align_summary, all, index_only };
	TASK task = TASK::all;
	BlastFormat blastFormat = BlastFormat::TABULAR;

	// methods
private:
	/*
	 * main method of this class.
	 * Parses the command line options, validates, and sets the class member variables
	 */
	void process(int argc, char**argv, bool dryrun);
	void validate();
	void validate_idxdir(); // called from validate
	void validate_kvdbdir(); // called from validate
	void validate_readb_dir(); // called from validate
	void validate_aligned_pfx();
	void validate_other_pfx();

	/*
	 * Resolve a user-supplied option name (any number of leading dashes, dashes
	 * and underscores interchangeable) to its entry in the 'options' table.
	 * Returns nullptr if the name is not a known option.
	 */
	const OptDesc* find_opt(const std::string& name) const;
	/*
	 * Rewrite argv into the form cxxopts accepts: canonical names, single dash
	 * for single-character names, double dash otherwise, and values attached
	 * with '='. Unknown dashed tokens are passed through unchanged so that
	 * cxxopts reports them as unrecognised.
	 */
	std::vector<std::string> normalize_argv(int argc, char** argv) const;
	/*
	 * Handle 'sortmerna <option> -h'. Returns true if a single-option help
	 * request was recognized and printed.
	 */
	bool print_single_option_help(int argc, char** argv);

	void opt_reads(const std::string& val);
	void opt_reads_gz(char **argv, int& narg);
	void opt_ref(const std::string& val);
	void opt_aligned(const std::string &val);
	void opt_other(const std::string &val);
	void opt_denovo_otu(const std::string &val);
	void opt_otu_map(const std::string &val);
	void opt_print_all_reads(const std::string &val);
	void opt_pid(const std::string &val);
	void opt_paired(const std::string& val);
	void opt_paired_in(const std::string &val);
	void opt_paired_out(const std::string &val);
	void opt_out2(const std::string& val);
	void opt_sout(const std::string& val);
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
	void opt_a(const std::string& val);
	void opt_e(const std::string& val); // opt_e_Evalue
	void opt_F(const std::string& val); // opt_F_ForwardOnly
	void opt_R(const std::string& val); // opt_R_ReverseOnly
	void opt_h(const std::string& val); // '-h' compact option list
	void opt_help(const std::string& val); // '--help' full reference
	void opt_v(const std::string& val); // opt_v_Verbose
	void opt_N(const std::string& val); // opt_N_MatchAmbiguous
	void opt_workdir(const std::string& path);
	void opt_kvdb(const std::string& path);
	void opt_idxdir(const std::string& path); // see help_idxdir
	void opt_readb(const std::string& path);
	void opt_dbg_level(const std::string& val);

	// ref interval m L max_pos v h  // indexing options
	void opt_interval(const std::string &val);
	void opt_m(const std::string &val);
	void opt_L(const std::string &val);
	void opt_max_pos(const std::string &val);
	void opt_readfeed(const std::string& val);
	void opt_score_split(const std::string& val);
	void opt_flush_delay(const std::string& val);
	/*
	 * true: 1,yes,Yes,Y,y,T,t, false: 0,No,NO,no,N,n,F,f
	*/
	void opt_zip_out(const std::string& val);
	void opt_index(const std::string& val); // help_index
	void opt_align(const std::string& val); // TODO: may be no need for this  20210207
	void opt_filter(const std::string& val); // TODO: may be no need for this  20210207

	void opt_dbg_put_db(const std::string& opt);
	void opt_max_read_len(const std::string& val);

	std::string to_string();
	std::string to_bin_string();
	//void store_to_db(KeyValueDatabase& kvdb);

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

	// container for options passed to the program. Keyed by canonical name;
	// iteration order is alphabetical, which is the option processing order.
	std::multimap<std::string, std::string> mopt;

	// OPTIONS table - specifies all possible options
	const std::array<OptDesc, 58> options = {{
		{OPT_REF,            "PATH",        COMMON,      true,  OptArg::REQUIRED, "",           help_ref,             help_ref_long,             &Runopts::opt_ref},
		{OPT_READS,          "PATH",        COMMON,      true,  OptArg::REQUIRED, "",           help_reads,           help_reads_long,           &Runopts::opt_reads},
		{OPT_WORKDIR,        "PATH",        COMMON,      false, OptArg::REQUIRED, "USRDIR/sortmerna/run/", help_workdir, help_workdir_long,      &Runopts::opt_workdir},
		{OPT_KVDB,           "PATH",        COMMON,      false, OptArg::REQUIRED, "WORKDIR/kvdb",  help_kvdb,         help_kvdb_long,            &Runopts::opt_kvdb},
		{OPT_IDXDIR,         "PATH",        COMMON,      false, OptArg::REQUIRED, "WORKDIR/idx",   help_idxdir,       "",                        &Runopts::opt_idxdir},
		{OPT_READB,          "PATH",        COMMON,      false, OptArg::REQUIRED, "WORKDIR/readb/",help_readb,        help_readb_long,           &Runopts::opt_readb},
		{OPT_READFEED,       "INT",         COMMON,      false, OptArg::REQUIRED, "0",          help_readfeed,        help_readfeed_long,        &Runopts::opt_readfeed},
		{OPT_FASTX,          "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_fastx,           "",                        &Runopts::opt_fastx},
		{OPT_SAM,            "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_sam,             "",                        &Runopts::opt_sam},
		{OPT_SQ,             "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_SQ,              "",                        &Runopts::opt_SQ},
		{OPT_BLAST,          "STR",         COMMON,      false, OptArg::REQUIRED, "1",          help_blast,           help_blast_long,           &Runopts::opt_blast},
		{OPT_ALIGNED,        "PATH",        COMMON,      false, OptArg::OPTIONAL, "WORKDIR/out/aligned", help_aligned, help_aligned_long,        &Runopts::opt_aligned},
		{OPT_OTHER,          "PATH",        COMMON,      false, OptArg::OPTIONAL, "WORKDIR/out/other",   help_other,   help_other_long,          &Runopts::opt_other},
		{OPT_NUM_ALIGNMENTS, "INT",         COMMON,      false, OptArg::REQUIRED, "1",          help_num_alignments,  help_num_alignments_long,  &Runopts::opt_num_alignments},
		{OPT_NO_BEST,        "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_no_best,         help_no_best_long,         &Runopts::opt_no_best},
		{OPT_MIN_LIS,        "INT",         COMMON,      false, OptArg::REQUIRED, "2",          help_min_lis,         help_min_lis_long,         &Runopts::opt_min_lis},
		{OPT_PRINT_ALL_READS,"BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_print_all_reads, "",                        &Runopts::opt_print_all_reads},
		{OPT_PAIRED,         "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_paired,          help_paired_long,          &Runopts::opt_paired},
		{OPT_PAIRED_IN,      "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_paired_in,       help_paired_in_long,       &Runopts::opt_paired_in},
		{OPT_PAIRED_OUT,     "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_paired_out,      help_paired_out_long,      &Runopts::opt_paired_out},
		{OPT_OUT2,           "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_out2,            help_out2_long,            &Runopts::opt_out2},
		{OPT_SOUT,           "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_sout,            help_sout_long,            &Runopts::opt_sout},
		{OPT_ZIP_OUT,        "STR/BOOL",    COMMON,      false, OptArg::OPTIONAL, "-1",         help_zip_out,         help_zip_out_long,         &Runopts::opt_zip_out},
		{OPT_MATCH,          "INT",         COMMON,      false, OptArg::REQUIRED, "2",          help_match,           "",                        &Runopts::opt_match},
		{OPT_MISMATCH,       "INT",         COMMON,      false, OptArg::REQUIRED, "-3",         help_mismatch,        "",                        &Runopts::opt_mismatch},
		{OPT_GAP_OPEN,       "INT",         COMMON,      false, OptArg::REQUIRED, "5",          help_gap_open,        "",                        &Runopts::opt_gap_open},
		{OPT_GAP_EXT,        "INT",         COMMON,      false, OptArg::REQUIRED, "2",          help_gap_ext,         "",                        &Runopts::opt_gap_ext},
		{OPT_E,              "DOUBLE",      COMMON,      false, OptArg::REQUIRED, "1e-5",       help_e,               help_e_long,               &Runopts::opt_e},
		{OPT_F,              "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_F,               "",                        &Runopts::opt_F},
		{OPT_N,              "INT",         COMMON,      false, OptArg::REQUIRED, "0",          help_N,               "",                        &Runopts::opt_N},
		{OPT_R,              "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_R,               "",                        &Runopts::opt_R},
		{OPT_SCORE_SPLIT,    "BOOL",        COMMON,      false, OptArg::NONE,     "False",      help_score_split,     help_score_split_long,     &Runopts::opt_score_split},
		{OPT_MAX_READ_LEN,   "INT",         COMMON,      false, OptArg::OPTIONAL, std::to_string(MAX_READ_LEN), help_max_read_len, "",           &Runopts::opt_max_read_len},
		{OPT_ID,             "DOUBLE",      OTU_PICKING, false, OptArg::REQUIRED, "0.97",       help_id,              "",                        &Runopts::opt_id},
		{OPT_COVERAGE,       "DOUBLE",      OTU_PICKING, false, OptArg::REQUIRED, "0.97",       help_coverage,        "",                        &Runopts::opt_coverage},
		{OPT_DENOVO_OTU,     "BOOL",        OTU_PICKING, false, OptArg::NONE,     "False",      help_denovo_otu,      help_denovo_otu_long,      &Runopts::opt_denovo_otu},
		{OPT_OTU_MAP,        "BOOL",        OTU_PICKING, false, OptArg::NONE,     "False",      help_otu_map,         help_otu_map_long,         &Runopts::opt_otu_map},
		{OPT_PASSES,         "INT,INT,INT", ADVANCED,    false, OptArg::REQUIRED, "L,L/2,3",    help_passes,          "",                        &Runopts::opt_passes},
		{OPT_EDGES,          "INT",         ADVANCED,    false, OptArg::REQUIRED, "4",          help_edges,           help_edges_long,           &Runopts::opt_edges},
		{OPT_NUM_SEEDS,      "INT",         ADVANCED,    false, OptArg::REQUIRED, "2",          help_num_seeds,       "",                        &Runopts::opt_num_seeds},
		{OPT_FULL_SEARCH,    "BOOL",        ADVANCED,    false, OptArg::NONE,     "False",      help_full_search,     help_full_search_long,     &Runopts::opt_full_search},
		{OPT_PID,            "BOOL",        ADVANCED,    false, OptArg::NONE,     "False",      help_pid,             "",                        &Runopts::opt_pid},
		{OPT_A,              "INT",         ADVANCED,    false, OptArg::REQUIRED, "numCores",   help_a,               help_a_long,               &Runopts::opt_a},
		{OPT_THREADS,        "INT",         ADVANCED,    false, OptArg::REQUIRED, "2",          help_threads,         "",                        &Runopts::opt_threads},
		{OPT_INDEX,          "INT",         INDEXING,    false, OptArg::OPTIONAL, "1",          help_index,           help_index_long,           &Runopts::opt_index},
		{OPT_L,              "DOUBLE",      INDEXING,    false, OptArg::REQUIRED, "18",         help_L,               "",                        &Runopts::opt_L},
		{OPT_M,              "DOUBLE",      INDEXING,    false, OptArg::REQUIRED, "3072",       help_m,               "",                        &Runopts::opt_m},
		{OPT_V,              "BOOL",        INDEXING,    false, OptArg::NONE,     "True",       help_v,               "",                        &Runopts::opt_v},
		{OPT_INTERVAL,       "INT",         INDEXING,    false, OptArg::REQUIRED, "1",          help_interval,        "",                        &Runopts::opt_interval},
		{OPT_MAX_POS,        "INT",         INDEXING,    false, OptArg::REQUIRED, "10000",      help_max_pos,         help_max_pos_long,         &Runopts::opt_max_pos},
		{OPT_H,              "BOOL",        HELP,        false, OptArg::NONE,     "",           help_h,               "",                        &Runopts::opt_h},
		{OPT_HELP,           "BOOL",        HELP,        false, OptArg::NONE,     "",           help_help,            "",                        &Runopts::opt_help},
		{OPT_VERSION,        "BOOL",        HELP,        false, OptArg::NONE,     "",           help_version,         "",                        &Runopts::opt_version},
		{OPT_DBG_PUT_DB,     "BOOL",        DEVELOPER,   false, OptArg::NONE,     "False",      help_dbg_put_db,      "",                        &Runopts::opt_dbg_put_db},
		{OPT_CMD,            "BOOL",        DEVELOPER,   false, OptArg::NONE,     "False",      help_cmd,             "",                        &Runopts::opt_cmd},
		{OPT_TASK,           "INT",         DEVELOPER,   false, OptArg::REQUIRED, "4",          help_task,            help_task_long,            &Runopts::opt_task},
		{OPT_DBG_LEVEL,      "INT",         DEVELOPER,   false, OptArg::REQUIRED, "0",          help_dbg_level,       help_dbg_level_long,       &Runopts::opt_dbg_level},
		{OPT_FLUSH_DELAY,    "INT",         DEVELOPER,   false, OptArg::REQUIRED, "60",         help_flush_delay,     help_flush_delay_long,     &Runopts::opt_flush_delay}
	}};
	// ~options table
}; // ~struct Runopts
// ~options.cpp

} // namespace sortmerna
