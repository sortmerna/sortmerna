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

/*
 * FILE: options.hpp
 * Created: Aug 19, 2017 Sat
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

/*! @brief Maximum length of input reads
    (not limited to this length algorithmically)
*/
const uint64_t MAX_READ_LEN = 30000;

// global constants
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
OPT_LOG = "log",
OPT_NUM_ALIGNMENTS = "num_alignments",
OPT_NO_BEST = "no-best",
OPT_MIN_LIS = "min_lis",
OPT_PRINT_ALL_READS = "print_all_reads",
OPT_PAIRED = "paired",
OPT_PAIRED_IN = "paired_in",
OPT_PAIRED_OUT = "paired_out",
OPT_OUT2 = "out2",
OPT_SOUT = "sout",
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
OPT_MAX_POS = "max_pos",
OPT_READS_FEED = "reads_feed",  // TODO: on hold
OPT_ZIP_OUT = "zip-out",
OPT_INDEX = "index",
OPT_ALIGN = "align",  // TODO: on hold
OPT_FILTER = "filter",  // TODO: on hold
OPT_DBG_LEVEL = "dbg-level",
OPT_MAX_READ_LEN = "max_read_len";

// help strings
const std::string \
help_header =
"  Usage:   sortmerna -ref FILE [-ref FILE] -reads FWD_READS [-reads REV_READS] [OPTIONS]:\n"
"  -------------------------------------------------------------------------------------------------------------\n"
"  | option            type-format           description                                          default      |\n"
"  -------------------------------------------------------------------------------------------------------------\n",
help_ref = 
	"Reference file (FASTA) absolute or relative path.\n\n"
	"       Use mutliple times, once per a reference file\n\n",
help_reads = 
	"Raw reads file (FASTA/FASTQ/FASTA.GZ/FASTQ.GZ).\n\n"
	"       Use twice for files with paired reads.\n"
	"       The file extensions are Not important. The program automatically\n"
	"       recognizes the file format as flat/compressed, fasta/fastq\n\n",
help_aligned = 
	"Aligned reads file prefix [dir/][pfx]       WORKDIR/out/aligned\n\n"
	"       Directory and file prefix for aligned output i.e. each\n"
	"       output file goes into the specified directory with the given prefix.\n"
	"       The appropriate extension: (fasta|fastq|blast|sam|etc) is automatically added.\n"
	"       Both 'dir' and 'pfx' are optional.\n"
	"       The 'dir' can be a relative or an absolute path.\n"
	"       If 'dir' is not specified, the output is created in the WORKDIR/out/\n"
	"       If 'pfx' is not specified, the prefix 'aligned' is used\n"
	"       Examples:\n"
	"       '-aligned $MYDIR/dir_1/dir_2/1' -> $MYDIR/dir_1/dir_2/1.fasta\n"
	"       '-aligned dir_1/apfx'           -> $PWD/dir_1/apfx.fasta\n"
	"       '-aligned dir_1/'               -> $PWD/aligned.fasta\n"
	"       '-aligned apfx'                 -> $PWD/apfx.fasta\n"
	"       '-aligned  (no argument)'       -> WORKDIR/out/aligned.fasta\n\n",
help_other = 
	"Non-aligned reads file prefix [dir/][pfx]   WORKDIR/out/other\n\n"
	"       Directory and file prefix for non-aligned output i.e. each\n"
	"       output file goes into the specified directory with the given prefix.\n"
	"       The appropriate extension: (fasta|fastq|blast|sam|etc) is automatically added.\n"
	"       Must be used with '" + OPT_FASTX + "'.\n"
	"       Both 'dir' and 'pfx' are optional.\n"
	"       The 'dir' can be a relative or an absolute path.\n"
	"       If 'dir' is not specified, the output is created in the WORKDIR/out/\n"
	"       If 'pfx' is not specified, the prefix 'other' is used\n"
	"       Examples:\n"
	"       '-other $MYDIR/dir_1/dir_2/1' -> $MYDIR/dir_1/dir_2/1.fasta\n"
	"       '-other dir_1/apfx'           -> $PWD/dir_1/apfx.fasta\n"
	"       '-other dir_1/'               -> $PWD/dir_1/other.fasta\n"
	"       '-other apfx'                 -> $PWD/apfx.fasta\n"
	"       '-other  (no argument)'       -> aligned_out/other.fasta\n"
	"                                        i.e. the same output directory\n"
	"                                        as used for aligned output\n\n",
help_fastx = 
	"Output aligned reads into FASTA/FASTQ file",
help_workdir = 
	"Workspace directory                         USRDIR/sortmerna/run/\n\n"
	"       Default structure: WORKDIR/\n"
	"                              idx/   (References index)\n"
	"                              kvdb/  (Key-value storage for alignments)\n"
	"                              out/   (processing output)\n"
	"                              readb/ (pre-processed reads/index)\n\n",
help_kvdb =
	"Directory for Key-value database            WORKDIR/kvdb\n\n"
	"       KVDB is used for storing the alignment results.\n\n",
help_idxdir =
	"Directory for storing Reference index.      WORKDIR/idx\n\n",
help_readb = 
	"Storage for pre-processed reads             WORKDIR/readb/\n\n"
	"       Directory storing the split reads, or the random access index of compressed reads\n\n",
	//"       Use with '" + OPT_READS_FEED + "'\n\n",
help_sam = 
	"Output SAM alignment for aligned reads.\n\n",
help_SQ = 
	"Add SQ tags to the SAM file\n\n",
help_blast = 
	"output alignments in various Blast-like formats\n\n"
	"       Sample values: '0'                    - pairwise\n"
	"                      '1'                    - tabular (Blast - m 8 format)\n"
	"                      '1 cigar'              - tabular + column for CIGAR\n"
	"                      '1 cigar qcov'         - tabular + columns for CIGAR and query coverage\n"
	"                      '1 cigar qcov qstrand' - tabular + columns for CIGAR, query coverage,\n"
	"                                               and strand\n\n",
help_dbg_put_db = 
	"",
help_log = 
	"Output overall statistics.                              True\n"
	"                                            TODO: remove\n",
help_num_alignments = 
	"Positive integer (INT >=0).\n\n"
	"       If used with '-" + OPT_NO_BEST + "' reports first INT alignments per read reaching\n"
	"       E-value threshold, which allows to lower the CPU time and memory use.\n"
	"       Otherwise outputs INT best alignments.\n"
	"       If INT = 0, all alignments are output\n\n",

help_no_best = 
	"Disable best alignments search                          False\n\n"
	"       The 'best' alignment is the highest scoring alignment out of All alignments of a read,\n"
	"       and the read can potentially be aligned (reaching E-value threshold) to multiple reference\n"
	"       sequences.\n"
	"       By default the program searches for best alignments i.e. performs an exhaustive search\n"
	"       over all references. Using '-" + OPT_NO_BEST + "' will make the program to search just\n"
	"       the first N alignments, where N is set using '-"+ OPT_NUM_ALIGNMENTS + "' i.e. 1 by default.\n\n",

help_min_lis = 
	"Search only alignments that have the LIS                2\n"
	"                                            of at least N seeds long\n\n"
	"       LIS stands for Longest Increasing Subsequence. It is computed using seeds, which\n"
	"       are k-mers common to the read and the reference sequence. Sorted sequences of such seeds\n"
	"       are used to filter the candidate references prior performing the Smith-Waterman alignment.\n\n",

help_print_all_reads = 
	"Output null alignment strings for non-aligned reads     False\n"
	"                                            to SAM and/or BLAST tabular files\n",
help_paired =
	"Flags paired reads                                      False\n\n"
	"        If a single reads file is provided, use this option to indicate\n"
	"        the file contains interleaved paired reads when neither\n"
	"        '" + OPT_PAIRED_IN + "' | '" + OPT_PAIRED_OUT + "' | '"+ OPT_OUT2 + "' | '" + OPT_SOUT + "' are specified.\n\n",
help_paired_in = 
	"Flags the paired-end reads as Aligned,                  False\n"
	"                                            when either of them is Aligned.\n\n"
	"        With this option both reads are output into Aligned FASTA/Q file\n"
	"        Must be used with '" + OPT_FASTX + "'.\n"
	"        Mutually exclusive with '" + OPT_PAIRED_OUT + "'.\n\n",

help_paired_out = 
	"Flags the paired-end reads as Non-aligned,              False\n"
	"                                            when either of them is non-aligned.\n\n"
	"        With this option both reads are output into Non-Aligned FASTA/Q file\n"
	"        Must be used with '" + OPT_FASTX + "'.\n"
	"        Mutually exclusive with '" + OPT_PAIRED_IN + "'.\n\n",

help_out2 =
	"Output paired reads into separate files.                False\n\n"
	"       Must be used with '" + OPT_FASTX + "'.\n"
	"       If a single reads file is provided, this options implies interleaved paired reads\n"
	"       When used with '"+ OPT_SOUT + "', four (4) output files for aligned reads will be generated:\n"
	"       'aligned-paired-fwd, aligned-paired-rev, aligned-singleton-fwd, aligned-singleton-rev'.\n"
	"       If '" + OPT_OTHER + "' option is also used, eight (8) output files will be generated.\n\n",

help_sout =
	"Separate paired and singleton aligned reads.            False\n\n"
	"       To be used with '" + OPT_FASTX + "'.\n"
	"       If a single reads file is provided, this options implies interleaved paired reads\n"
	"       Cannot be used with '" + OPT_PAIRED_IN + "' | '" + OPT_PAIRED_OUT + "'\n\n",

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
	"E-value threshold.                                      1\n\n"
	"       Defines the 'statistical significance' of a local alignment.\n"
	"       Exponentially correllates with the Minimal Alignment score.\n"
	"       Higher E-values (100, 1000, ...) cause More reads to Pass the alignment threshold\n\n",

help_v = 
	"Produce verbose output when building the index          True\n",

help_id = 
	"%%id similarity threshold (the alignment                0.97\n"
	"                                            must still pass the E-value threshold).\n",

help_coverage = 
	"%%query coverage threshold (the alignment must          0.97\n"
	"                                            still pass the E-value threshold)\n",

help_denovo_otu = 
	"Output FASTA file with 'de novo' reads                  False\n\n"
	"       Read is 'de novo' if its alignment score passes E-value threshold, but both the identity\n"
	"       '-" + OPT_ID + "', and the '-" + OPT_COVERAGE + "' are below their corresponding thresholds\n"
	"       i.e. ID < %%id and COV < %%cov\n\n",

help_otu_map = 
	"Output OTU map (input to QIIME's make_otu_table.py).    False\n"
	"                                            Cannot be used with '" + OPT_NO_BEST + " because\n"
	"                                            the grouping is done around the best alignment'\n",
help_passes = 
	"Three intervals at which to place the seed on           L,L/2,3\n"
	"                                             the read (L is the seed length)\n",
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
	"Processing Task                                         4\n\n"
	"       Possible values: 0 - align. Only perform alignment\n"
	"                        1 - post-processing (log writing)\n"
	"                        2 - generate reports\n"
	"                        3 - align and post-process\n"
	"                        4 - all\n\n",
help_a = 
	"DEPRECATED in favour of '-threads'. Number of           numCores\n"
	"                                            processing threads to use.\n"
	"                                            Automatically redirects to '-threads'\n",
help_threads = 
	"Number of Processing threads to use                     2\n",
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
	"Indexing: the amount of memory (in Mbytes) for          3072\n"
	"                                            building the index.\n",

help_L = 
	"Indexing: seed length.                                  18\n",

help_max_pos = 
	"Indexing: maximum (integer) number of positions to      1000\n"
	"                                            store for each unique L-mer.\n"
	"                                            If 0 - all positions are stored.\n",

//help_reads_feed = 
//	"Method of accessing the reads by the                    0\n"
//	"                                            reads processors\n\n"
//	"       0 - Split reads. Reads files are split into parts equal the number of processing threads\n"
//	"       1 - FUTURE: Lockless queue. Reads are put into a lockless queue\n"
//	"                   to be popped by the processing threads\n"
//	"       3 - FUTURE: Random access to the compresssed reads files\n"
//	"       4 - FUTURE: combination of the random access and the lockless queue\n\n",

help_zip_out =
	"Controls the output compression                        '-1'\n\n"
	"       By default the report files are produced in the same format as the input i.e.\n"
	"       if the reads files are compressed (gz), the output is also compressed.\n"
	"       The default behaviour can be overriden by using '-" + OPT_ZIP_OUT + "'.\n"
	"       The possible values: '1/true/t/yes/y'\n"
	"                            '0/false/f/no/n'\n"
	"                            '-1' (the same format as input - default)\n"
	"       The values are Not case sensitive i.e. 'Yes, YES, yEs, Y, y' are all OK\n"
	"       Examples:\n"
	"       '-" + OPT_READS + " freads.gz -" + OPT_ZIP_OUT + " n' : generate flat output when the input is compressed\n"
	"       '-" + OPT_READS + " freads.flat -" + OPT_ZIP_OUT + "' : compress the output when the input files are flat\n\n",

help_index =
    "Build reference database index                          2\n\n"
	"       By default when this option is not used, the program checks the reference index and\n"
	"       builds it if not already existing.\n"
	"       This can be changed by using '-" + OPT_INDEX + "' as follows:\n"
	"       '-" + OPT_INDEX + " 0' - skip indexing. If the index does not exist, the program will terminate\n"
	"                                and warn to build the index prior performing the alignment\n"
	"       '-" + OPT_INDEX + " 1' - only perform the indexing and terminate\n"
	"       '-" + OPT_INDEX + " 2' - the default behaviour, the same as when not using this option at all\n\n",

help_dbg_level =
	"Debug level                                             0\n\n"
	"      Controls verbosity of the execution trace. Default value of 0 corresponds to\n"
	"      the least verbose output.\n"
	"      The highest value currently is 2.\n\n",

help_max_read_len =
	"Maximum allowed read length                             " + std::to_string(MAX_READ_LEN) + "\n\n"

//help_align =
//    "Perform the alignment                                   False\n\n"
//	"       Search a single best alignment per read\n\n",
//
//help_filter =
//    "Perform the filtering                                   False\n\n"
//	"       Search for a single first found alignment per read\n\n"
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
	/*
	 * verify the string is a valid option
	*/
	bool is_option(const std::string & opt);

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
	bool is_paired_out = false; // '--paired_out' both paired-end reads go in 'other' fasta/q file. Only Fasta/q and De-novo reporting.
	bool is_out2 = false; // 20200127 output paired reads into separate files. Issue 202
	bool is_sout = false; // 20210105 separate singletons and paired
	bool is_otu_map = false; // OPT_OTU_MAP was selected i.e. output OTU map (input to QIIME's make_otu_table.py)
	bool is_denovo = false; // output file with reads matching database < %%id (set using --id) and < %%cov (set using --coverage)
	bool is_log = true; // OPT_LOG was selected i.e. output overall statistics. TODO: remove this option, always generate.
	bool is_print_all_reads = false; // '--print_all_reads' output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files
	bool is_sam = false; // OPT_SAM was specified. output SAM alignment (for aligned reads only)
	bool is_SQ = false; // OPT_SQ add SQ tags to the SAM file
	bool is_blast = false; // OPT_BLAST was specified
	bool is_fastx = false; // OPT_FASTX was selected i.e. output FASTA/FASTQ file (for aligned and/or rejected reads)
	bool is_other = false; // OPT_OTHER was selected i.e. flags to produce 'other' file
	bool is_verbose; // OPT_V was selected (indexing)
	bool is_pid = false; // add pid to output file names
	bool is_cmd = false; // start interactive session
	bool is_dbg_put_kvdb = false; // if True - do Not put records into Key-value DB. Debugging Memory Consumption.
	int  findex = 2; // 0 (don't build index) | 1 (only build index) | 2 (default - build index if not present)
	bool is_align = false;
	bool is_filter = false;

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

	int num_read_thread = 1;     // number of threads reading the Reads file.
	int num_write_thread = 1;    // number of threads writing to Key-value database
	int num_proc_thread = 2;     // number of threads to use for alignment
	int num_read_thread_pp = 1;  // number of post-processing read threads
	int num_proc_thread_pp = 1;  // number of post-processing processor threads
	int num_read_thread_rep = 1; // number of report reader threads
	int num_proc_thread_rep = 1; // number of report processor threads
	int dbg_level = 0; // lowest debug level - minimal info.

	int queue_size_max = 1000; // max number of Reads in the Read and Write queues. 10 works OK.
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
	long gap_open = 5; // '--gap_open' SW penalty (positive integer) for introducing a gap
	long gap_extension = 2; // '--gap_ext' SW penalty (positive integer) for extending a gap
	int score_N = 0; // '-N' SW penalty for ambiguous letters (N's)
	FEED_TYPE feed_type = FEED_TYPE::SPLIT_READS; // OPT_READS_FEED

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

	const std::string dbkey = "run_options";
	const std::string IDX_DIR  = "idx";
	const std::string KVDB_DIR = "kvdb";
	const std::string OUT_DIR  = "out";
	const std::string READB_DIR = "readb";

	enum ALIGN_REPORT { align, summary, report, alnsum, all, index_only };
	ALIGN_REPORT alirep = ALIGN_REPORT::all;
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
	void opt_idxdir(const std::string& path); // see help_idxdir
	void opt_readb(const std::string& path);
	void opt_dbg_level(const std::string& val);

	// ref tmpdir interval m L max_pos v h  // indexing options
	void opt_tmpdir(const std::string &val);
	void opt_interval(const std::string &val);
	void opt_m(const std::string &val);
	void opt_L(const std::string &val);
	void opt_max_pos(const std::string &val);
	void opt_reads_feed(const std::string& val);
	/*
	 * true: 1,yes,Yes,Y,y,T,t, false: 0,No,NO,no,N,n,F,f
	*/
	void opt_zip_out(const std::string& val);
	void opt_index(const std::string& val); // help_index
	void opt_align(const std::string& val); // TODO: may be no need for this  20210207
	void opt_filter(const std::string& val); // TODO: may be no need for this  20210207

	void opt_default(const std::string& opt);
	void opt_dbg_put_db(const std::string& opt);
	void opt_unknown(char** argv, int& narg, char* opt);
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

	// container for options passed to the program
	std::multimap<std::string, std::string> mopt;

	// OPTIONS Map - specifies all possible options
	const std::array<opt_6_tuple, 54> options = {
		std::make_tuple(OPT_REF,            "PATH",        COMMON,      true,  help_ref, &Runopts::opt_ref),
		std::make_tuple(OPT_READS,          "PATH",        COMMON,      true,  help_reads, &Runopts::opt_reads),
		//std::make_tuple(OPT_ALIGN,          "BOOL",        COMMON,      true,  help_align, &Runopts::opt_align),
		//std::make_tuple(OPT_FILTER,         "BOOL",        COMMON,      true,  help_filter, &Runopts::opt_filter),
		std::make_tuple(OPT_WORKDIR,        "PATH",        COMMON,      false, help_workdir, &Runopts::opt_workdir),
		std::make_tuple(OPT_KVDB,           "PATH",        COMMON,      false, help_kvdb, &Runopts::opt_kvdb),
		std::make_tuple(OPT_IDXDIR,         "PATH",        COMMON,      false, help_idxdir, &Runopts::opt_idxdir),
		std::make_tuple(OPT_READB,          "PATH",        COMMON,      false, help_readb, &Runopts::opt_readb),
		std::make_tuple(OPT_FASTX,          "BOOL",        COMMON,      false, help_fastx, &Runopts::opt_fastx),
		std::make_tuple(OPT_SAM,            "BOOL",        COMMON,      false, help_sam, &Runopts::opt_sam),
		std::make_tuple(OPT_SQ,             "BOOL",        COMMON,      false, help_SQ, &Runopts::opt_SQ),
		std::make_tuple(OPT_BLAST,          "STR",         COMMON,      false, help_blast, &Runopts::opt_blast),
		std::make_tuple(OPT_ALIGNED,        "STR/BOOL",    COMMON,      false, help_aligned, &Runopts::opt_aligned),
		std::make_tuple(OPT_OTHER,          "STR/BOOL",    COMMON,      false, help_other, &Runopts::opt_other),
		std::make_tuple(OPT_NUM_ALIGNMENTS, "INT",         COMMON,      false, help_num_alignments, &Runopts::opt_num_alignments),
		std::make_tuple(OPT_NO_BEST,        "BOOL",        COMMON,      false, help_no_best, &Runopts::opt_no_best),
		std::make_tuple(OPT_MIN_LIS,        "INT",         COMMON,      false, help_min_lis, &Runopts::opt_min_lis),
		std::make_tuple(OPT_PRINT_ALL_READS,"BOOL",        COMMON,      false, help_print_all_reads, &Runopts::opt_print_all_reads),
		std::make_tuple(OPT_PAIRED,         "BOOL",        COMMON,      false, help_paired, &Runopts::opt_paired),
		std::make_tuple(OPT_PAIRED_IN,      "BOOL",        COMMON,      false, help_paired_in, &Runopts::opt_paired_in),
		std::make_tuple(OPT_PAIRED_OUT,     "BOOL",        COMMON,      false, help_paired_out, &Runopts::opt_paired_out),
		std::make_tuple(OPT_OUT2,           "BOOL",        COMMON,      false, help_out2, &Runopts::opt_out2),
		std::make_tuple(OPT_SOUT,           "BOOL",        COMMON,      false, help_sout, &Runopts::opt_sout),
		std::make_tuple(OPT_ZIP_OUT,        "STR/BOOL",    COMMON,      false, help_zip_out, &Runopts::opt_zip_out),
		std::make_tuple(OPT_MATCH,          "INT",         COMMON,      false, help_match, &Runopts::opt_match),
		std::make_tuple(OPT_MISMATCH,       "INT",         COMMON,      false, help_mismatch, &Runopts::opt_mismatch),
		std::make_tuple(OPT_GAP_OPEN,       "INT",         COMMON,      false, help_gap_open, &Runopts::opt_gap_open),
		std::make_tuple(OPT_GAP_EXT,        "INT",         COMMON,      false, help_gap_ext, &Runopts::opt_gap_ext),
		std::make_tuple(OPT_E,              "DOUBLE",      COMMON,      false, help_e, &Runopts::opt_e),
		std::make_tuple(OPT_F,              "BOOL",        COMMON,      false, help_F, &Runopts::opt_F),
		std::make_tuple(OPT_N,              "BOOL",        COMMON,      false, help_N, &Runopts::opt_N),
		std::make_tuple(OPT_R,              "BOOL",        COMMON,      false, help_R, &Runopts::opt_R),
		std::make_tuple(OPT_MAX_READ_LEN,   "INT",         COMMON,      false, help_max_read_len, &Runopts::opt_max_read_len),
		//std::make_tuple(OPT_READS_FEED,     "INT",         COMMON,      false, help_reads_feed, &Runopts::opt_reads_feed),
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
		std::make_tuple(OPT_INDEX,          "INT",         INDEXING,    false, help_index, &Runopts::opt_index),
		std::make_tuple(OPT_L,              "DOUBLE",      INDEXING,    false, help_L, &Runopts::opt_L),
		std::make_tuple(OPT_M,              "DOUBLE",      INDEXING,    false, help_m, &Runopts::opt_m),
		std::make_tuple(OPT_V,              "BOOL",        INDEXING,    false, help_v, &Runopts::opt_v),
		std::make_tuple(OPT_INTERVAL,       "INT",         INDEXING,    false, help_interval, &Runopts::opt_interval),
		std::make_tuple(OPT_MAX_POS,        "INT",         INDEXING,    false, help_max_pos, &Runopts::opt_max_pos),
		std::make_tuple(OPT_H,              "BOOL",        HELP,        false, help_h, &Runopts::opt_h),
		std::make_tuple(OPT_VERSION,        "BOOL",        HELP,        false, help_version, &Runopts::opt_version),
		std::make_tuple(OPT_DBG_PUT_DB,     "BOOL",        DEVELOPER,   false, help_dbg_put_db, &Runopts::opt_dbg_put_db),
		std::make_tuple(OPT_CMD,            "BOOL",        DEVELOPER,   false, help_cmd, &Runopts::opt_cmd),
		std::make_tuple(OPT_TASK,           "INT",         DEVELOPER,   false, help_task, &Runopts::opt_task),
		std::make_tuple(OPT_DBG_LEVEL,      "INT",         DEVELOPER,   false, help_dbg_level, &Runopts::opt_dbg_level)
		//std::make_tuple(OPT_THREP,          "INT:INT",     DEVELOPER,   false, help_threp, &Runopts::opt_threp)
	};
	// ~map options
}; // ~struct Runopts
// ~options.cpp
