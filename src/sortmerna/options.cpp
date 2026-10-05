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
 * FILE: options.cpp
 * Created: Jun 07, 2018 Thu
 */

 // TODO: BUG: if SMR headers moved down after 3rd party, 'timeval' struct gets 'redefined' - compiler error. That's a header mess bug.

 // 3rd party
#include "zlib.h"
#include "cxxopts.hpp"

#include "version.h"
#include "build_version.h"
#include "options.hpp"
#include "common.hpp"
#include "izlib.hpp"
#include "kvdb.hpp"

 // standard
#include <limits>
#include <dirent.h>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring> // strerror, strrchr, memcpy, strcpy, strpbrk
#include <fcntl.h>
#include <functional> // std::invoke
#include <filesystem>
#include <algorithm>
#include <iomanip>


#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

// forward

namespace sortmerna {

void about();
//void help();
std::string get_user_home(); // util.cpp
std::string get_basename(const std::string& file); // util.cpp
std::streampos filesize(const std::string& file); // util.cpp

Runopts::Runopts(int argc, char** argv, bool dryrun)
{
	process(argc, argv, dryrun);
}

/*
 * Guard for the handlers that convert the value with stoi/stod and have no
 * meaningful behaviour without one. An option name given with no value yields
 * an empty string (see normalize_argv), so without this the conversion would
 * throw std::invalid_argument with no indication of which option is at fault.
 */
static void require_value(const std::string& opt, const std::string& val, const std::string& example)
{
	if (val.empty()) {
		SMR_THROW("Option '--", opt, "' requires a value e.g. '--", opt, " ", example, "'");
	}
}

/* 
 * validates file (be it plain or zipped) and sets corresponding options
 * For paired reads files use --reads file1 --reads file2
 */
void Runopts::opt_reads(const std::string &file)
{
	std::stringstream ss;

	auto numread = mopt.count(OPT_READS);
	auto readcnt = readfiles.size();

	INFO("Processing reads file [", readcnt + 1, "] out of total [", numread, "] files");

	if (file.size() == 0)
	{
		SMR_THROW(help_reads);
	}

	// check file exists
	auto fpath = std::filesystem::path(file);
	auto fpath_a = std::filesystem::path(); // absolute path

	if (std::filesystem::exists(fpath))
	{
		fpath_a = fpath;
	}
	else if (fpath.is_relative())
	{
		INFO("File  [", file, "] appears to be a relative path. Trying to locate in current directory and in the working directory ...");
		fpath_a = std::filesystem::current_path() / file; // std::filesystem::path(file) makes no difference
		if (std::filesystem::exists(fpath_a))
		{
			std::cout << "File [" << fpath_a << "] exists" << std::endl;
		}
		else
		{
			INFO("File  [", file, "] has not been found in the current directory. Trying the working directory ...");
			if (!std::filesystem::exists(workdir / file))
			{
				INFO("Could not locate File [", file, "] neither at current path [", std::filesystem::current_path(), "] nor in Workdir [", workdir, "]");
				SMR_THROW(help_reads);
			}
		}
	}
	else
	{
		SMR_THROW("The file [", file, "] is not an existing/valid absolute or relative path\n", help_reads);
	}

	// check the file can be read
	std::ifstream ifs(fpath_a, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open())
	{
		SMR_THROW("Failed to open file [" , fpath_a , "]");
	}
	else {
		// verify the file is not empty (unless only building the index, task = index_only)
		// NOTE: opt_task may not have run yet, so read the raw task value from mopt.
		auto task_it = mopt.find(OPT_TASK);
		bool is_index_only = (task_it != mopt.end())
			&& (std::stoi(task_it->second) == static_cast<int>(TASK::index_only));
		ifs.seekg(0, std::ios_base::end);
		auto fsize = ifs.tellg();
		if (fsize <= 0 && !is_index_only)
		{
			SMR_THROW("The reads file [", fpath_a, "] is empty");
		}
		ifs.close();
		have_reads = true;
		readfiles.push_back(fpath_a.generic_string());
	}
} // ~Runopts::opt_reads

void Runopts::opt_ref(const std::string &refpath)
{
	auto numref = mopt.count(OPT_REF);
	auto refcnt = indexfiles.size();

	INFO("Processing reference [", refcnt + 1, "] out of total [", numref, "] references");

	if (refpath.size() == 0)
	{
		SMR_THROW(help_ref);
	}

	// check file exists and can be read
	auto fpath = std::filesystem::path(refpath);
	auto fpath_a = std::filesystem::path(); // absolute path

	if (std::filesystem::exists(fpath))
	{
		fpath_a = fpath;
	}
	else if (fpath.is_relative())
	{
		INFO("File  [" , refpath , "] appears to be a relative path. Trying to locate in current directory and in the working directory ...");
		fpath_a = std::filesystem::current_path() / refpath; // std::filesystem::path(file) makes no difference
		if (std::filesystem::exists(fpath_a))
		{
			INFO("File [" , fpath_a , "] exists");
		}
		else
		{
			INFO("File  [" , refpath , "] has not been found in the current directory. Trying the working directory ...");
			if (!std::filesystem::exists(workdir / refpath))
			{
				SMR_THROW("Could not locate File [" , refpath , "] neither at current path [" ,
					std::filesystem::current_path()	, "] nor in Workdir [" , workdir , "]\n", help_ref);
			}
		}
	}
	else
	{
		SMR_THROW("The file " , refpath , " is not an existing/valid absolute or relative path\n", help_ref);
	}

	// check files are readable
	std::ifstream ifstr(fpath_a);
	if (!ifstr.is_open() || !ifstr.good()) {
		SMR_THROW("Cannot read file [" , refpath , "]");
	}
	else {
		INFO("File " , std::filesystem::absolute(refpath) , " exists and is readable");
	}

	if (ifstr.is_open())
		ifstr.close();

	// check index file names are distinct
	for (int i = 0; i < (int)indexfiles.size(); i++)
	{
		if ((indexfiles[i].first).compare(refpath) == 0)
		{
			WARN("Reference file (" , refpath , ") has been entered more than once. Ignoring redundant enties");
		}
	}
	// add reference file absolute path. The index file will be set during index initialization.
	indexfiles.push_back(std::pair<std::string, std::string>(fpath_a.generic_string(), ""));
} // ~Runopts::opt_ref

/* 
 * possible values:
 *   no arg          WORKDIR/out/aligned
 *   .               PWD/aligned
 *   pfx             PWD/out/pfx
 *   dir1/pfx        PWD/dir1/pfx   use WORKDIR instead of PWD? - this would be a non-standard behaviour - confusing.
 *   ./dir1/pfx      PWD/dir1/pfx
 *   /dir1/dir2/pfx  /dir1/dir2/pfx
 *   /dir1/dir2/     /dir1/dir2/aligned
 */
void Runopts::opt_aligned(const std::string &file)
{
	// -aligned specified without argument
	if (file.size() == 0)
	{
		INFO("Directory and Prefix for the aligned output was not provided. Using default dir/pfx: 'WORKDIR/out/aligned'");
	}
	else {
		std::filesystem::path fpath = file;
		if (!fpath.empty()) {
			if (fpath.has_filename()) {
				aligned_pfx = fpath; // prefix is non-empty - use it
			}
			else {
				aligned_pfx = fpath / OPT_ALIGNED; // prefix doesn't specify the file name e.g. 'dir_1/' not 'dir_1/pfx_1'
			}
		}
	}
} // ~Runopts::opt_aligned

/**
 * depends on '--fastx'
 */
void Runopts::opt_other(const std::string &file)
{
	auto cnt = mopt.count(OPT_FASTX);
	if (cnt == 0)
	{
		SMR_THROW("Option '" + OPT_OTHER + "' can only be used together with '"+ OPT_FASTX + "' option.");
	}

	if (file.size() == 0)
	{
		std::cout << STAMP << OPT_OTHER << " was specified without argument. Will use default Directory and Prefix for the non-aligned output." << std::endl;
	}
	else {
		std::filesystem::path fpath = file;
		if (!fpath.empty()) {
			if (fpath.has_filename()) {
				other_pfx = fpath; // prefix is non-empty - use it
			}
			else {
				other_pfx = fpath / OPT_OTHER;
			}
		}
	}
	is_other = true;
} // ~Runopts::opt_other

void Runopts::opt_denovo_otu(const std::string& val)
{
	is_denovo = true;
} // ~Runopts::opt_de_novo_otu

void Runopts::opt_otu_map(const std::string &val)
{
	is_otu_map = true;
} // ~Runopts::opt_otu_map

void Runopts::opt_print_all_reads(const std::string &val)
{
	is_print_all_reads = true;
} // ~Runopts::optPrintAllReads

void Runopts::opt_pid(const std::string &val)
{
	is_pid = true;
} // ~Runopts::optPid

void Runopts::opt_paired(const std::string& val)
{
	std::stringstream ss;
	auto numread = mopt.count(OPT_READS);
	if (numread > 1) {
		ss << STAMP << "'" << OPT_PAIRED << "' " 
			"can only be used with a single reads file to indicate it holds paired reads.\n"
			"However option '" << OPT_READS << "' was specified [" << numread << "] times";
		SMR_THROW(ss.str());
	}
	is_paired = true;
} // ~Runopts::optPaired

void Runopts::opt_paired_in(const std::string &val)
{
	is_paired_in = true;
} // ~Runopts::optPairedIn

void Runopts::opt_paired_out(const std::string &val)
{
	is_paired_out = true;
} // ~Runopts::optPairedOut

/**
 * two 'reads' options has to be provided
 *
 */
void Runopts::opt_out2(const std::string& val)
{
	is_out2 = true;
} // ~Runopts::opt_out2

void Runopts::opt_sout(const std::string& val)
{
	is_sout = true;
}

void Runopts::opt_match(const std::string &val)
{
	std::stringstream ss;
	if (val.size() == 0)
	{
		ss << STAMP << "'" << OPT_MATCH << "' " << "requires a positive integer as input e.g. 2";
		SMR_THROW(ss.str());
	}
	// set match
	if (!match_set)
	{
		match = atoi(val.data());
		match_set = true;
	}
	else
	{
		ss.str("");
		ss << STAMP << "'" << OPT_MATCH << "' " << "[INT] has been set twice, please verify your choice";
		SMR_THROW(ss.str());
	}
} // ~Runopts::optMatch

void Runopts::opt_mismatch(const std::string &val)
{
	if (val.size() == 0)
	{
		SMR_THROW("--mismatch [INT] requires a negative integer input (ex. --mismatch -2)");
	}

	// set mismatch
	if (!mismatch_set)
	{
		mismatch = std::stoi(val);
		if (mismatch > 0)
		{
			SMR_THROW("--mismatch [INT] takes a negative integer (ex. --mismatch -2)");
		}
		mismatch_set = true;
	}
	else
	{
		print_help();
		SMR_THROW("--mismatch [INT] has been set twice, please verify your choice");
	}
} // ~Runopts::opt_mismatch

void Runopts::opt_gap_open(const std::string &val)
{
	if (val.size() == 0)
	{
		SMR_THROW("--gap_open [INT] requires a positive integer as input (ex. --gap_open 5)");
	}

	// set gap open
	if (!gap_open_set)
	{
		gap_open = atoi(val.data());
		if (gap_open < 0)
		{
			SMR_THROW("--gap_open [INT] requires a positive integer as input (ex. --gap_open 5)");
		}
		gap_open_set = true;
	}
	else
	{
		print_help();
		SMR_THROW("--gap_open [INT] has been set twice, please verify your choice");
	}
} // ~Runopts::opt_gap_open

void Runopts::opt_gap_ext(const std::string &val)
{
	if (val.size() == 0)
	{
		SMR_THROW("--gap_ext [INT] requires a positive integer as input (ex. --gap_ext 2)");
	}
	// set gap extend
	if (!gap_ext_set)
	{
		gap_extension = atoi(val.data());
		if (gap_extension < 0)
		{
			SMR_THROW("--gap_ext [INT] requires a positive integer as input (ex. --gap_ext 2)");
		}
		gap_ext_set = true;
	}
	else
	{
		print_help();
		SMR_THROW("--gap_ext [INT] has been set twice, please verify your choice");
	}
} // ~Runopts::opt_gap_ext

void Runopts::opt_num_seeds(const std::string& val)
{
	if (val.size() == 0)
	{
		SMR_THROW("--num_seeds [INT] requires a positive integer as input (ex. --num_seeds 6)");
	}

	char* end = 0;
	num_seeds = (int)strtol(val.data(), &end, 10); // convert to integer
	if (num_seeds <= 0)
	{
		SMR_THROW("--num_seeds [INT] requires a positive integer (>0) as input (ex. --num_seeds 6)");
	}

} // ~Runopts::opt_num_seeds

/* --fastx */
void Runopts::opt_fastx(const std::string &val)
{
	if (is_fastx)
	{
		SMR_THROW("--fastx has already been set once.");
	}
	else
	{
		is_fastx = true;
	}
} // ~Runopts::opt_fastx

void Runopts::opt_sam(const std::string &val)
{
	if (is_sam)
	{
		WARN("'" , OPT_SAM , "' has already been set. Ignoring second flag.");
	}
	else
	{
		is_sam = true;
	}
} // ~Runopts::opt_sam

void Runopts::opt_blast(const std::string &val)
{
	std::stringstream ss;

	if (is_blast)
	{
		ss << STAMP << "'" << OPT_BLAST << "' [STRING] has already been set once. Ignoring the second value";
		WARN(ss.str());
	}

	// split blast options into vector by space
	std::istringstream iss(val);
	do
	{
		std::string s;
		iss >> s;
		blastops.push_back(s);
	} while (iss);

	// remove the end of file entry
	blastops.pop_back();
	bool blast_human_readable = false;
	std::vector<std::string> supported_opts;
	supported_opts.push_back("0");
	supported_opts.push_back("1");
	supported_opts.push_back("cigar");
	supported_opts.push_back("qstrand");
	supported_opts.push_back("qcov");
	// check user options are supported
	for (uint32_t i = 0; i < blastops.size(); i++)
	{
		bool match_found = false;
		std::string opt = blastops[i];
		std::vector<std::string>::iterator it;
		for (it = supported_opts.begin(); it != supported_opts.end(); ++it)
		{
			if (opt.compare(*it) == 0)
			{
				if (opt.compare("0") == 0) {
					blast_human_readable = true;
					blastFormat = BlastFormat::REGULAR;
				}
				else if (opt.compare("1") == 0) {
					blastFormat = BlastFormat::TABULAR;
				}
				match_found = true;
				break;
			}
		}
		if (!match_found)
		{
			SMR_THROW("'" , opt , "' is not supported in --blast [STRING].");
		}
	}
	// more than 1 field with blast human-readable format given
	if (blast_human_readable && (blastops.size() > 1))
	{
		SMR_THROW("for human-readable format, --blast [STRING] can only contain a single field '0'.");
	}
	// both human-readable and tabular format options have been chosen
	if (blast_human_readable && blastFormat == BlastFormat::TABULAR)
	{
		SMR_THROW("--blast [STRING] can only have one of the options '0' (human-readable) or '1' (tabular).");
	}

	is_blast = true;
} // ~Runopts::opt_blast

void Runopts::opt_min_lis(const std::string &val)
{
	std::stringstream ss;
	if (val.size() == 0)
	{
		SMR_THROW("'", OPT_MIN_LIS, "' [INT] requires a positive integer as input e.g. 2. (if 0 - all high scoring reference sequences are searched)");
	}

	// min_lis_gv has already been set
	if (is_min_lis)
	{
		SMR_THROW("'" , OPT_MIN_LIS , "' [INT] has been set twice, please verify your choice.");
	}
	else
	{
		if ((sscanf(val.data(), "%d", &min_lis) != 1) || min_lis < 0)
		{
			SMR_THROW("'", OPT_MIN_LIS, "' [INT] requires a positive integer as input e.g. 2. If 0, all high scoring reference sequences are searched)");
		}
		is_min_lis = true;
	}
} // ~Runopts::opt_min_lis

void Runopts::opt_no_best(const std::string &val)
{
	is_best = false;
	INFO("'", OPT_NO_BEST, "' was selected. Disabling Best search.");
} // ~Runopts::opt_best

void Runopts::opt_num_alignments(const std::string &val)
{
	if (val.size() == 0)
	{
		SMR_THROW("'", OPT_NUM_ALIGNMENTS, "' [INT] requires a posistive integer as input e.g. 2. If 0, all alignments are output.");
	}

	if (is_num_alignments)
	{
		SMR_THROW("'", OPT_NUM_ALIGNMENTS, "' [INT] has been set twice, please verify your parameters.");
	}

	// set number of alignments to output reaching the E-value
	auto ii = std::stoi(val);
	if (ii < 0)
	{
		SMR_THROW("'", OPT_NUM_ALIGNMENTS, "' requires a positive integer as input e.g. 2. If 0, all alignments are output.");
	}
	else {
		num_alignments = ii;
	}

	is_num_alignments = true;
} // ~Runopts::opt_num_alignments

void Runopts::opt_edges(const std::string &val)
{
	require_value(OPT_EDGES, val, "4");
	// --edges is already set
	if (edges_set)
	{
		SMR_THROW(" : --edges [INT]%% has already been set once.");
	}

	char *end = 0;
	// find if % sign exists
	if (val.find_first_of("%") != std::string::npos)
		is_as_percent = true;

	// convert to integer
	edges = std::stoi(val);

	if (edges < 1 || edges > 10)
	{
		SMR_THROW(" : --edges [INT]%% requires a positive integer between 0-10 as input (ex. --edges 4).");
	}
} // ~Runopts::optEdges

void Runopts::opt_full_search(const std::string &val)
{
	if (full_search_set)
	{
		WARN("Options '--full_search' has been set more than once. Only the last flag is considered.");
	}
	full_search_set = true;
	is_full_search = true;
} // ~Runopts::opt_full_search

void Runopts::opt_SQ(const std::string &val)
{
	if (is_SQ)
	{
		SMR_THROW(" : BOOL --SQ has been set twice, please verify your choice.");
	}

	is_SQ = true;
} // ~Runopts::optSQ

void Runopts::opt_passes(const std::string &val)
{
	require_value(OPT_PASSES, val, "18,9,3");
    auto countpass = mopt.count(OPT_PASSES);
	if (countpass > 1)
	{
		SMR_THROW("'", OPT_PASSES, "' [INT,INT,INT] has been set more than once");
	}

    INFO("'", OPT_PASSES, "' option has been set. Using the provided skip lengths.");
    // Parse comma-separated tokens "INT,INT,INT". Only parse and validate here;
    // the values are stashed in 'skiplengths_user' and replicated per-index at
    // the end of process(), once 'indexfiles' is populated. Options are processed
    // in alphabetical order ('passes' before 'ref'), so 'indexfiles' is still
    // empty at this point.
    std::vector<uint32_t> skiplenvec;
    std::size_t count = 0;
    for (std::size_t pos_from = 0; pos_from <= val.size(); )
    {
    	auto pos = val.find(",", pos_from);
    	auto tok = val.substr(pos_from, pos == std::string::npos ? std::string::npos : pos - pos_from);
    	pos_from = (pos == std::string::npos) ? val.size() + 1 : pos + 1;
    	if (++count > 3)
    	{
    		SMR_THROW("Exactly 3 integers has to be provided with '" , OPT_PASSES , "' [INT,INT,INT]");
    	}
    	auto skiplen = std::stoi(tok);
    	if (skiplen > 0)
    		skiplenvec.emplace_back(skiplen);
    	else
    	{
    		SMR_THROW("All three integers in '", OPT_PASSES, "' [INT,INT,INT] must contain positive integers where 0 < INT < (shortest read length).");
    	}
    } // ~for
    if (skiplenvec.size() != 3)
    {
        SMR_THROW("Exactly 3 integers has to be provided with '" , OPT_PASSES , "' [INT,INT,INT]");
    }
    skiplengths_user = skiplenvec;
	passes_set = true;
} // ~Runopts::opt_passes

void Runopts::opt_id(const std::string &val)
{
	// % id
	if (min_id < 0)
	{
		if ((sscanf(val.data(), "%lf", &min_id) != 1) ||
			(min_id < 0) || (min_id > 1))
		{
			SMR_THROW(" : --id [DOUBLE] must be a positive float with value 0<=id<=1.");
		}
	}
	else
	{
		SMR_THROW(" : --id [DOUBLE] has been set twice, please verify your command parameters.");
	}
} // ~Runopts::opt_id

void Runopts::opt_coverage(const std::string &val)
{
	// % query coverage
	if (min_cov < 0)
	{
		if ((sscanf(val.data(), "%lf", &min_cov) != 1) ||
			(min_cov < 0) || (min_cov > 1))
		{
			SMR_THROW(" : --coverage [DOUBLE] must be a positive float with value 0<=id<=1.");
		}
	}
	else
	{
		SMR_THROW(" : --coverage [DOUBLE] has been set twice, please verify your command parameters.");
	}
} // ~Runopts::opt_coverage

void Runopts::opt_version(const std::string &val)
{
	std::cout << std::endl
		<< "SortMeRNA version " << SORTMERNA_MAJOR << "." << SORTMERNA_MINOR << "." << SORTMERNA_PATCH << std::endl
		<< "Build Date: " << sortmerna_build_compile_date << std::endl
		<< sortmerna_build_git_sha << std::endl
		<< sortmerna_build_git_date << std::endl;
	throw smr_exit_requested("--version");
} // ~Runopts::opt_version

void Runopts::opt_e(const std::string &val)
{
	std::stringstream ss;

	// E-value
	if (val.size() == 0)
	{
		SMR_THROW(": -e [DOUBLE] requires a positive double as input (ex. --e 1e-5)");
	}

	if (evalue < 0)
	{
		sscanf(val.data(), "%lf", &evalue);  // default 1e-5 is set in validate()
		if (evalue < 0)
		{
			SMR_THROW(" : -e [DOUBLE] requires a positive double as input (ex. --e 1e-5)");
		}
	}
	else
	{
		SMR_THROW(" : -e [DOUBLE] has been set twice, please verify your command parameters.");
	}
} // ~Runopts::opt_e

void Runopts::opt_F(const std::string &val)
{
	// only forward strand
	if (!is_forward)
	{
		is_forward = true;
	}
	else
	{
		SMR_THROW(" : BOOL -F has been set more than once, please check your command parameters.");
	}
} // ~Runopts::opt_F_ForwardOnly

void Runopts::opt_R(const std::string &val)
{
	// only reverse strand
	if (!is_reverse)
	{
		is_reverse = true;
	}
	else
	{
		SMR_THROW(" : BOOL '-R' has been set more than once, please check your command parameters.");
	}
} // ~Runopts::opt_R

void Runopts::opt_h(const std::string &val)
{
	about();
	print_help();
	throw smr_exit_requested("-h");
} // ~Runopts::opt_h

void Runopts::opt_help(const std::string &val)
{
	about();
	print_help_full();
	throw smr_exit_requested("--help");
} // ~Runopts::opt_help

void Runopts::opt_v(const std::string &val)
{
	is_verbose = true;
} // ~Runopts::opt_v

void Runopts::opt_N(const std::string &val)
{
	require_value(OPT_N, val, "0");
	if (!match_ambiguous_N)
	{
		match_ambiguous_N = true;
		score_N = std::stoi(val);
	}
	else
	{
		SMR_THROW(": BOOL -N has been set more than once, please check your command parameters.");
	}
} // ~Runopts::opt_N

/* Number Processor threads to use */
void Runopts::opt_a(const std::string &val)
{
	WARN("Option '-a' was Deprecated in favour of '-threads'. Redirecting to 'threads' ...");
	//if (val.size() == 0)
	//{
	//	ERR(": -a [INT] requires an integer for number of Processor threads (ex. -a 8)");
	//	exit(EXIT_FAILURE);
	//}
	opt_threads(val);
	//num_proc_thread = std::stoi(val);
} // ~Runopts::opt_a_numProcThreads

/* Number of threads to use 
 * @param val INT[:INT[:INT]] e.g. 8 | 8:1 | 8:1:1 i.e. takes at least one integer value
 *                                 |     |       |_ write threads (future)
 *                                 |     |_ read threads (future)
 *                                 |_ processor threads
 */
void Runopts::opt_threads(const std::string &val)
{
	std::string msg = "'-threads INT' requires an integer for number of "
		"Processing threads e.g. '-threads 8'. Default value equals the Number of CPUs or CPU cores";
	std::istringstream strm(val);
	std::string tok;

	if (val.size() == 0)
	{
		SMR_THROW(msg);
	}
	//else
	//{
	//	auto n = std::count(val.begin(), val.end(), ':');
	//	if (n != 2) {
	//		ERR(msg);
	//		exit(EXIT_FAILURE);
	//	}
	//}

	for (int i = 0; std::getline(strm, tok, ':'); ++i)
	{
		switch (i)
		{
		case 0: num_proc_thread = std::stoi(tok); break;
		case 1: 
			WARN("Using more than a single Read thread is reserved for future implementation. Now using 1 thread");
			//num_read_thread = std::stoi(tok); 
			break;
		case 2:
			WARN("Using more than a single Write thread is reserved for future implementation. Now using 1 thread");
			//num_write_thread = std::stoi(tok); 
			break;
		}
	}
} // ~Runopts::opt_threads


void Runopts::opt_dbg_put_db(const std::string& val)
{
	is_dbg_put_kvdb = true;
}

  /* Processing task */
void Runopts::opt_task(const std::string &val)
{
	require_value(OPT_TASK, val, "4");
	int task_num = std::stoi(val);

	if (task_num > 5)
	{
		SMR_THROW("Option '", OPT_TASK, "' can only take values in range [0..5] Provided value is [", task_num , "'");
	}

	switch (task_num)
	{
	case 0: task = align; break;
	case 1: task = summary; break;
	case 2: task = report; break;
	case 3: task = align_summary; break;
	case 4: task = all; break;
	case 5: task = index_only; break;
	}
} // ~Runopts::opt_task

void Runopts::opt_dbg_level(const std::string& val)
{
	require_value(OPT_DBG_LEVEL, val, "1");
	dbg_level = std::stoi(val);
	if (dbg_level > 2) {
		INFO("provided value: ", dbg_level, " for '", OPT_DBG_LEVEL, "' is too high. Using the currently highest execution trace verbosity of 2");
		dbg_level = 2;
	}
}

void Runopts::opt_flush_delay(const std::string& val)
{
	require_value(OPT_FLUSH_DELAY, val, "60");
	int v = std::stoi(val);
	if (v <= 0) {
		SMR_THROW("'", OPT_FLUSH_DELAY, "' must be a positive integer (seconds); got: ", val);
	}
	flush_delay = static_cast<unsigned>(v);
}

// interactive session '--cmd'
void Runopts::opt_cmd(const std::string &val)
{
	is_cmd = true;
} // ~Runopts::optInteractive

void Runopts::opt_workdir(const std::string &path)
{
	std::stringstream ss;
	if (path.size() == 0)
	{
		SMR_THROW("'" , OPT_WORKDIR, "' option takes an argument - a directory path. None was provided.\n", help_workdir);
	}
	else {
		workdir = path;
		INFO("Using WORKDIR: ", std::filesystem::absolute(workdir), " as specified");
	}
}

void Runopts::opt_kvdb(const std::string& path) {
	if (path.size() == 0)
	{
		SMR_THROW("'" , OPT_KVDB, "' option takes an argument - a directory path. None was provided.\n" , help_kvdb);
	}
	else {
		kvdbdir = path;
		INFO("Using KVDB dir: [" , std::filesystem::absolute(kvdbdir) , " as specified");
	}
}

void Runopts::opt_idxdir(const std::string& path) {
	if (path.size() == 0)
	{
		SMR_THROW("'" , OPT_IDXDIR, "' option takes an argument - a directory path. None was provided.\n" , help_kvdb);
	}
	else {
		idxdir = path;
		INFO("Using IDX dir: [" , std::filesystem::absolute(idxdir) , " as specified");
	}
}

void Runopts::opt_readb(const std::string& path) {
	if (path.size() == 0)
	{
		SMR_THROW("'", OPT_READB, "' option takes an argument - a directory path. None was provided.\n", help_readb);
	}
	else {
		readb_dir = path;
		INFO("Using Reads DB dir: [", std::filesystem::absolute(readb_dir), " as specified");
	}
}

void Runopts::opt_interval(const std::string &val)
{
	auto count = mopt.count(OPT_INTERVAL);
	if (count > 1)
	{
		WARN("Option 'interval' entered [" , count , "] times. Only the last value will be used. " , help_interval);
	}

	if (val.size() == 0)
	{
		WARN("Option 'interval' given without value. Using default: ", interval);
	}
	else
	{
		interval = std::stoi(val);
	}
}

/**
 * Max size of an index file/part 
 */
void Runopts::opt_m(const std::string &val)
{
	auto count = mopt.count(OPT_M);
	if (count > 1)
	{
		WARN(" Option 'm' entered [" , count , "] times. Only the last value will be used. " , help_m);
	}

	if (val.size() == 0)
	{
		WARN("Option 'interval' given without value. Using default: 3072 MB");
	}
	else
	{
		max_file_size = std::stod(val);
	}
}

void Runopts::opt_L(const std::string &val)
{
	auto count = mopt.count(OPT_L);
	if (count > 1)
	{
		SMR_THROW(" Option '", OPT_L, "' entered [", count, "] times. Only a single one is expected.\n", 
                "\tHelp: ", help_L);
	}

	if (val.size() == 0)
	{
		WARN("Option 'L' given without value. Using default: " , seed_win_len);
	}
	else
	{
		int lnwin_t = std::stoi(val);

		if (lnwin_t <= 0 || lnwin_t % 2 == 1 || lnwin_t < 8 || lnwin_t > 26)
		{
			WARN("Option '", OPT_L, "' takes a Positive Even integer between 8 and 26 inclusive"
				" e.g. 10, 12, 14, .. , 20. Provided value: ", lnwin_t, " Default will be used: ", seed_win_len);
		}
		else
		{
			seed_win_len = lnwin_t;
		}
	}
} // ~Runopts::opt_L

void Runopts::opt_max_pos(const std::string &val)
{
	auto count = mopt.count(OPT_MAX_POS);
	if (count > 1)
	{
		WARN("Option '" , OPT_MAX_POS , "' entered [" , count , "] times. Only the last value will be used\n"
			, "\tHelp: " , help_max_pos);
	}

	if (val.size() == 0)
	{
		WARN("Options 'max_pos' takes a positive integer e.g. 250. Using default: " , max_pos);
	}
	else
	{
		max_pos = std::stoi(val);
	}
}

void Runopts::opt_readfeed(const std::string& val)
{
	require_value(OPT_READFEED, val, "0");
	FEED_TYPE ftype = static_cast<FEED_TYPE>(std::stoi(val));

	if (ftype > FEED_TYPE::MAX)
	{
		SMR_THROW("Option '", OPT_READFEED, "' can only take values in range [0..", static_cast<int>(FEED_TYPE::MAX), "] Provided value is ['", val, "'");
	}

	feed_type = ftype;
} // ~Runopts::opt_readfeed

void Runopts::opt_zip_out(const std::string& val)
{
	const std::array<std::string, 5> yesvals = {"1", "y", "yes", "t", "true"};
	const std::array<std::string, 5> novals = {"0", "n", "no", "f", "false"};
	if (val.size() > 0) {
		if (val != "-1") {
			// to lowercase
			std::string valc(val);
			std::transform(valc.begin(), valc.end(), valc.begin(),
							[](unsigned char c) -> unsigned char { return std::tolower(c); });
			auto pval = std::find(std::begin(yesvals), std::end(yesvals), valc);
			if (pval != std::end(yesvals)) {
				zip_out = 1;
			}
			else {
				auto pval = std::find(std::begin(novals), std::end(novals), valc);
				if (pval != std::end(novals)) {
					zip_out = 0;
				}
			}

			if (zip_out != 0 && zip_out != 1) {
				WARN("'", OPT_ZIP_OUT, "' was provided with an unrecognized value: ", val, " Using default: ", zip_out);
			}
			else {
				INFO("using '", OPT_ZIP_OUT, "' with specified value ", val, "(", zip_out, ")");
			}
		}
		else {
			INFO("using '", OPT_ZIP_OUT, "' with default value ", val);
		}
	}
}

void Runopts::opt_index(const std::string& val)
{
	if (val.size() > 0) {
		auto ii = std::stoi(val);
		if (!(ii == 0 || ii == 1)) {
			WARN("'", OPT_INDEX, "' can only take integer values: [", 0, ", ", 1, "]. Provided value: ", val, " Using default: ", findex);
		}
		else {
			INFO("using '", OPT_INDEX, "' with specified value ", ii);
			findex = ii;
		}
	}
	else {
		INFO("No value was provided with '", OPT_INDEX, "'. Using default: ", findex, " If no index exists it will be built");
	}
}

void Runopts::opt_align(const std::string& val)
{
	is_align = true;
}

void Runopts::opt_filter(const std::string& val)
{
	is_filter = true;
}

void Runopts::opt_max_read_len(const std::string& val)
{
	if (val.size() > 0) {
		max_read_len = std::stoull(val);
		INFO("using '", OPT_MAX_READ_LEN, "' with specified value ", max_read_len);
	}
	else {
		INFO("using '", OPT_MAX_READ_LEN, "' with default value ", max_read_len);
	}
}

void Runopts::opt_score_split(const std::string& val)
{
	is_score_split = true;
}

/* 
 * called from validate
 */
void Runopts::validate_idxdir() {
	if (idxdir.empty()) {
		if (workdir.empty()) {
			INFO("'" , OPT_WORKDIR	, "' option was not provided. Using USERDIR as the location for the Index");
			workdir = std::filesystem::path(get_user_home()) / WORKDIR_DEF_SFX;
		}
		idxdir = workdir / IDX_DIR; // default
	}
	INFO("Using index directory: " , std::filesystem::absolute(idxdir));
	if (!std::filesystem::exists(idxdir)) {
		bool is_dir_ok = std::filesystem::create_directory(idxdir);
		if (!is_dir_ok) {
			SMR_THROW("Failed creating IDX directory: [" , std::filesystem::absolute(idxdir) , "]");
		}
		else {
			INFO("Created index directory - OK");
		}
	}
	else {
		if (std::filesystem::is_empty(idxdir)) {
			INFO("IDX directory: " , std::filesystem::absolute(idxdir) , " exists and is empty");
		}
		else {
			INFO("IDX directory: " , std::filesystem::absolute(idxdir) , " exists and is not empty");
		}
	}
}

/* 
 * called from validate
 */
void Runopts::validate_kvdbdir()
{
	if (kvdbdir.empty()) {
		if (workdir.empty()) {
			INFO("'", OPT_WORKDIR, "' option was not provided. Using USERDIR to set the working directory: ", workdir);
			workdir = std::filesystem::path(get_user_home()) / WORKDIR_DEF_SFX;
		}
		kvdbdir = workdir / KVDB_DIR; // default
	}

	INFO("Key-value DB location " , std::filesystem::absolute(kvdbdir));

	if (std::filesystem::exists(kvdbdir))
	{
		// dir exists and is empty
		if (std::filesystem::is_empty(kvdbdir))
		{
			// dir exists and empty -> use
			if (TASK::summary == task || TASK::report == task)
			{
				INFO("KVDB directory: " , std::filesystem::absolute(kvdbdir) , " is empty. OK to use.");
			}
		}
		else // not empty
		{
			// Non-empty kvdb is the resume path. restart::probe_or_init will
			// verify the opts/inputs fingerprints when the DB is opened and
			// either resume (matching) or abort with a diagnostic (mismatched).
			std::stringstream ss;
			for (auto& subpath : std::filesystem::directory_iterator(kvdbdir))
				ss << '\t' << subpath.path().filename() << '\n';
			INFO("KVDB directory ", std::filesystem::absolute(kvdbdir),
			     " is not empty. Will attempt auto-resume; existing contents:\n", ss.str());
		}
	}
	else
	{
		// dir does not exist -> try creating
		INFO("Creating KVDB directory: " , std::filesystem::absolute(kvdbdir));
		bool is_dir_ok = std::filesystem::create_directories(kvdbdir);
		if (!is_dir_ok) {
			SMR_THROW("Failed creating KVDB directory: " , std::filesystem::absolute(kvdbdir));
		}
	}
} // ~validate_kvdbdir

/*
  validate split reads directory
*/
void Runopts::validate_readb_dir() {
	const auto SR = "split reads";
	const auto SR_DIR = "split reads directory";
	if (readb_dir.empty()) {
		if (workdir.empty()) {
			INFO("'", OPT_WORKDIR, "' option was not provided. Using USERDIR as the location for the ",SR);
			workdir = std::filesystem::path(get_user_home()) / WORKDIR_DEF_SFX;
		}
		readb_dir = workdir / READB_DIR; // default
	}
	INFO("Using ",SR_DIR," : ", std::filesystem::absolute(readb_dir));
	if (!std::filesystem::exists(readb_dir)) {
		bool is_dir_ok = std::filesystem::create_directory(readb_dir);
		if (!is_dir_ok) {
			SMR_THROW("Failed creating ",SR_DIR," : ", std::filesystem::absolute(readb_dir));
		}
		else {
			INFO("Created ",SR_DIR," - OK");
		}
	}
	else {
		if (std::filesystem::is_empty(readb_dir)) {
			INFO(SR_DIR," : ", std::filesystem::absolute(readb_dir), " exists and is empty");
		}
		else {
			INFO(SR_DIR," : ", std::filesystem::absolute(readb_dir), " exists and is not empty");
		}
	}
} // ~Runopts::validate_readb_dir

void Runopts::validate_aligned_pfx() {
	if (aligned_pfx.empty()) {
		aligned_pfx = workdir / OUT_DIR / OPT_ALIGNED; // default output file
	}

	if (aligned_pfx.has_parent_path()) {
		// create directory if not exists
		if (!std::filesystem::exists(aligned_pfx.parent_path())) {
			INFO("Checking output directory: " , std::filesystem::absolute(aligned_pfx.parent_path()));
			bool is_dir_ok = std::filesystem::create_directories(aligned_pfx.parent_path());
			if (!is_dir_ok) {
				SMR_THROW("Failed creating output directory: ", aligned_pfx.parent_path().string());
			}
		}
	}
}

void Runopts::validate_other_pfx() {
	if (other_pfx.empty()) {
		if (aligned_pfx.empty()) {
			other_pfx = workdir / OUT_DIR / OPT_OTHER; // default output file
		}
		else {
			other_pfx = aligned_pfx.parent_path() / OPT_OTHER;
		}
	}

	if (other_pfx.has_parent_path()) {
		// create directory if not exists
		if (!std::filesystem::exists(other_pfx.parent_path())) {
			INFO("Checking output directory: " , std::filesystem::absolute(other_pfx.parent_path()));
			bool is_dir_ok = std::filesystem::create_directories(other_pfx.parent_path());
			if (!is_dir_ok) {
				SMR_THROW("Failed creating output directory: " , other_pfx.parent_path().string());
			}
		}
	}
}

void Runopts::process(int argc, char**argv, bool dryrun)
{
#if defined(__APPLE__)
	int sz[2] = { CTL_HW, HW_MEMSIZE };
	u_int namelen = sizeof(sz) / sizeof(sz[0]);
	uint64_t size;
	size_t len = sizeof(size);
	if (sysctl(sz, namelen, &size, &len, NULL, 0) < 0)
	{
		SMR_THROW("sysctl (main.cpp)");
	}
#endif

	if (argc == 1)
	{
		about();
		SMR_THROW("Missing required command options");
	}

	// store the command line
	for (int i = 0; i < argc; i++)
	{
		cmdline.append(argv[i]);
		cmdline.append(" ");
	}

	// 'sortmerna <option> -h' - help on a single option. Handled before anything
	// else because the bare option name is not a command line option itself.
	if (print_single_option_help(argc, argv)) {
		throw smr_exit_requested("single option help");
	}

	// Rewrite argv into the form cxxopts accepts, then let cxxopts tokenize it.
	// Every option is registered with an implicit value so that cxxopts never
	// consumes the token following an option name - normalize_argv has already
	// attached the values with '='. That keeps values that look like options
	// (e.g. '--mismatch -3') intact.
	auto norm = normalize_argv(argc, argv);
	std::vector<const char*> norm_argv;
	norm_argv.reserve(norm.size());
	for (const auto& tok : norm)
		norm_argv.push_back(tok.c_str());

	// The description and type are registered too, even though the help is
	// rendered from the 'options' table rather than by cxxopts (it has no
	// per-option help API and cannot express the required/default columns).
	// That keeps 'cxx.help()' usable as a fallback.
	cxxopts::Options cxx(argv[0], "");
	cxx.set_width(100);
	for (const auto& opt : options) {
		cxx.add_options(opt_category_name_map.at(opt.cat))(
			opt.name, opt.help,
			cxxopts::value<std::vector<std::string>>()->implicit_value(""), opt.type);
	}

	try {
		auto parsed = cxx.parse(static_cast<int>(norm_argv.size()), norm_argv.data());
		// 'arguments()' yields only the options actually given, in command line
		// order, one entry per occurrence. Feeding them into 'mopt' re-sorts them
		// alphabetically, which is the option processing order the handlers expect.
		for (const auto& kv : parsed.arguments())
			mopt.emplace(kv.key(), kv.value());
		for (const auto& extra : parsed.unmatched()) {
			SMR_THROW("Unexpected argument: '", extra, "'. Note that e.g. '-", OPT_REF,
				"' or '-", OPT_READS, "' have to be used with Each file. See 'sortmerna -h'");
		}
	}
	catch (const cxxopts::exceptions::exception& e) {
		SMR_THROW(e.what(), ". See 'sortmerna -h'");
	}

	// is this only Help options?
	bool is_help_opt = mopt.count(OPT_H) + mopt.count(OPT_HELP) + mopt.count(OPT_VERSION) > 0;

	// Answer help and version before the processing trace starts, so that the
	// output is not preceded by diagnostics. All three handlers throw
	// smr_exit_requested i.e. nothing below this point runs.
	if (is_help_opt)
	{
		for (const auto& opt : mopt) {
			if (opt.first == OPT_H || opt.first == OPT_HELP || opt.first == OPT_VERSION) {
				std::invoke(find_opt(opt.first)->fn, this, opt.second);
			}
		}
	}

	INFO("=== Options processing starts ... ===");

	if (!is_help_opt)
	{
		// check required options were provided
		for (const auto& opt : options)
		{
			if (opt.required && mopt.count(opt.name) == 0)
			{
				// NOTE: option-processing (incl. opt_task) hasn't run yet, so the
				// 'task' member still holds its default. Read the raw value from mopt.
				auto task_it = mopt.find(OPT_TASK);
				bool is_index_only = (task_it != mopt.end())
					&& (std::stoi(task_it->second) == static_cast<int>(TASK::index_only));
				if (opt.name == OPT_READS && is_index_only) {
					// reads are not required when only building the index
					continue;
				}
				else if (opt.name == OPT_READS) {
					SMR_THROW("Missing required flag: ", opt.name,
						". Use '", OPT_TASK, "' option with value '", TASK::index_only,
						"' if you intend to build index without processing reads. See 'sortmerna -h' for details.");
				}
				else {
					SMR_THROW("Missing required flag: ", opt.name, ". See 'sortmerna -h' for details.");
				}
			}
		}

		// process WORKDIR first (if it was specified) as other options depend on it
		if (mopt.count(OPT_WORKDIR) > 0) {
			auto wd_it = mopt.find(OPT_WORKDIR);
			if (wd_it != mopt.end())
			{
				std::string wdir = wd_it->second; // get the value
				mopt.erase(wd_it); // remove to prevent repeated processing further
				opt_workdir(wdir);
			}
		}
	}

	// loop through the rest of options
	for (auto opt : mopt)
	{
		INFO("Processing option: " , opt.first , " with value: " , opt.second);

		// cxxopts rejects unknown names, so the lookup always succeeds
		const OptDesc* desc = find_opt(opt.first);
		if (desc == nullptr) {
			SMR_THROW("Option: '", opt.first, "' is not recognized");
		}
		std::invoke(desc->fn, this, opt.second);
	}

	INFO("=== Options processing done ===\n");
	INFO("Alignment type: [best:", is_best,
			" num_alignments:", num_alignments,
			" min_lis:", min_lis, " seeds:", num_seeds,"]");

	if (!is_help_opt)
	{
		// Build per-index skip lengths now that 'indexfiles' is populated (opt_ref
		// ran in the loop above). Done here rather than in opt_passes because options
		// are processed in alphabetical order, so 'passes' is handled before 'ref'.
		if (skiplengths.empty())
		{
			std::vector<uint32_t> def;
			if (passes_set) {
				def = skiplengths_user; // user-supplied '--passes' [INT,INT,INT]
				INFO("Using '", OPT_PASSES, "' skip lengths: ", def[0], ",", def[1], ",", def[2]);
			}
			else if (mopt.count(OPT_L) > 0) {
				auto len = std::stoi(mopt.find(OPT_L)->second);
				def = { static_cast<uint32_t>(len), static_cast<uint32_t>(len / 2), 3 };
				INFO("'", OPT_PASSES, "' not set. Using '", OPT_L, "'-derived skip lengths: ",
					def[0], ",", def[1], ",", def[2]);
			}
			else {
				// Leave zeros; Refstats::load fills index-derived defaults.
				def = { 0, 0, 0 };
			}
			for (std::size_t i = 0; i < indexfiles.size(); ++i)
				skiplengths.push_back(def);
		}

		validate();
		//about(); // if we are here, args are OK, welcome the user
	}
} // ~Runopts::process

/**
 * Validate the options and setup defaults
 */
void Runopts::validate()
{
	validate_kvdbdir();
	validate_idxdir();
	if (feed_type == FEED_TYPE::SPLIT_READS)
		validate_readb_dir();
	validate_aligned_pfx(); // there is always some output like log => validate
	if (is_other) {
		validate_other_pfx();
	}

	// No output format has been chosen
	if (!(is_fastx || is_blast || is_sam || is_otu_map || is_denovo))
	{
		is_blast = true;
		INFO("No output format has been chosen (fastx|sam|blast|otu_map). Using default '" , OPT_BLAST , "'");
	}

	if (is_paired_in && is_paired_out)
	{
		SMR_THROW("Options '" , OPT_PAIRED_IN , "' and '" , OPT_PAIRED_OUT , "' are mutually exclusive. Please choose one or the other");
	}

	if (!is_paired) {
		is_paired = readfiles.size() == 2 || is_paired_in || is_paired_out;
	}

	if (is_out2 && !is_paired) {
		WARN("Option '", OPT_OUT2, 
			"' is Ignored because it can only be used with paired reads."
			" The reads are considered paired if either 2 reads files are supplied, or '", 
			OPT_PAIRED_IN, "', or '", OPT_PAIRED_OUT, "' is specified");
		is_out2 = false;
	}

	if (is_sout && !is_paired) {
		WARN("Option '", OPT_SOUT, "' is Ignored because it can only be used with paired reads."
			" for the purpose of '", OPT_SOUT, "' the reads are considered paired if either",
			" 2 reads files are supplied, or '", OPT_PAIRED, "' is specified with a single reads file");
		is_out2 = false;
	}

	if (is_sout && (is_paired_in || is_paired_out)) {
		SMR_THROW("Option '", OPT_SOUT,"' cannot be used when either '", OPT_PAIRED_IN,"' or '", OPT_PAIRED_OUT,"' is specified.");
	}

	// Options --paired_in and --paired_out can only be used with FASTA/Q output
	if (!is_fastx && (is_paired_in || is_paired_out))
	{
		INFO("Options '" , OPT_PAIRED_IN , "' and '" , OPT_PAIRED_OUT, 
			"' must be accompanied by option '" , OPT_FASTX , "'. Setting to true.");
		is_fastx = true;
	}

	// An OTU map can only be constructed with the single best alignment per read
	if (is_otu_map && !is_best)
	{
		SMR_THROW("'-", OPT_OTU_MAP, "' cannot be set together with '-", OPT_NO_BEST, "'.\n"
			"\tThe best alignment is required for constructing an OTU map.");
	}

	if (is_num_alignments && !(is_blast || is_sam || is_fastx))
	{
		SMR_THROW("'" , OPT_NUM_ALIGNMENTS, 
			"' [INT] has been set but no output format has been chosen (--blast | --sam | --fastx). Using default '", 
			OPT_BLAST , "'");
	}

	// Check gap extend score < gap open score
	if (gap_extension > gap_open)
	{
		SMR_THROW("--gap_ext [INT] must be less than --gap_open [INT].");
	}

	// Option --print_all_reads can only be used with Blast-like tabular
	// and SAM formats (not pairwise)
	if (is_print_all_reads && is_blast && blastFormat != BlastFormat::TABULAR)
	{
		SMR_THROW("--print_all_reads [BOOL] can only be used with BLAST-like");
	}

	// Option --min_lis [INT] can only accompany --best [INT]
	if (is_min_lis && is_num_alignments)
	{
		SMR_THROW("'" , OPT_MIN_LIS , "' [INT] and '" , OPT_NUM_ALIGNMENTS , "' [INT] cannot be set together.\n");
	}

	// Option --mis_lis INT accompanies --best INT, cannot be set alone
	if (is_min_lis && !is_best)
	{
		SMR_THROW("--min_lis [INT] must be set together with --best [INT].");
	}

	// %id and %coverage can only be used with --otu_map
	if ((min_id > 0 || min_cov > 0) && !is_otu_map)
	{
		SMR_THROW("--id [INT] and --coverage [INT] can only be used together with --otu_map.\n"
			"\tThese two options are used for constructing the OTU map\n"
			"\tby filtering alignments passing the E-value threshold.");
	}

	// if neither strand was selected for search, search both
	if (!is_forward && !is_reverse)
	{
		is_forward = true;
		is_reverse = true;
	}

	// default E-value 
	if (evalue < 0.0) {
        INFO("No E-value threshold provided. Using default: ", 1e-5);
        evalue = 1e-5;
	}

	// SW alignment parameters
	if (!match_set) {
		match = 2;
		match_set = true;
	}

	if (!mismatch_set) {
		mismatch = -3;
		mismatch_set = true;
	}

	if (!gap_open_set) {
		gap_open = 5;
		gap_open_set = true;
	}

	if (!gap_ext_set) {
		gap_extension = 2;
		gap_ext_set = true;
	}

	if (!match_ambiguous_N) 
		score_N = mismatch;

	// default method for searching alignments
	if (!is_best && !is_num_alignments)
	{
		// TODO: looks arbitrary. Why the alignment contolling options would depend on the output?
		// FASTA/FASTQ output, stop searching for alignments after the first match
		if (is_fastx && !(is_blast || is_sam || is_otu_map || is_log || is_denovo))
			num_alignments = 1;
		// output single best alignment from best candidate hits
		else
		{
			min_lis = 2;
		}
	}

	// default minimum LIS used for setting the number of
	// alignments to search prior to outputting --best INT
	if (is_best && !is_min_lis) 
		min_lis = 2;

	// default number of seed hits before searching for candidate LIS
	if (num_seeds < 0)
		num_seeds = 2;

	// default number of nucleotides to add to each edge
	// of an alignment region before extension
	if (edges < 0) 
		edges = 4;

	// activate heuristic for stopping search (of 1-error matches) after
	// finding 0-error match
	if (!full_search_set) 
		is_full_search = false;

	// default %id to keep alignment
	if (min_id < 0)
	{
		// if OTU-map is chosen, set default similarity to 0.97
		if (is_otu_map) min_id = 0.97;
		else min_id = 0;
	}

	// default %query coverage to keep alignment
	if (min_cov < 0)
	{
		// if OTU-map is chosen, set default coverage to 0.97
		if (is_otu_map) min_cov = 0.97;
		else min_cov = 0;
	}
} // ~Runopts::validate

/* 
 * human readable representation of the options
 */
std::string Runopts::to_string()
{
	return "TODO";
}

/** 
 * encoded options for storing in the Key-value database 
 */
std::string Runopts::to_bin_string()
{
	return "TODO";
}

/*
 * Map a user-supplied option name to its canonical spelling.
 * Dashes and underscores are interchangeable; canonical form uses dashes.
 * Case is significant - '-N' and '-n' are not the same option.
 */
static std::string canonicalize_opt_name(const std::string& name)
{
	std::string canon(name);
	std::replace(canon.begin(), canon.end(), '_', '-');
	return canon;
}

const Runopts::OptDesc* Runopts::find_opt(const std::string& name) const
{
	auto canon = canonicalize_opt_name(name);
	for (const auto& opt : options) {
		if (opt.name == canon)
			return &opt;
	}
	return nullptr;
}

/*
 * Greedy word wrap. Returns the lines without any indentation.
 */
static std::vector<std::string> wrap_text(const std::string& text, std::size_t width)
{
	std::vector<std::string> lines;
	std::istringstream iss(text);
	std::string word;
	std::string line;
	while (iss >> word) {
		if (line.empty()) {
			line = word;
		}
		else if (line.size() + 1 + word.size() <= width) {
			line += " " + word;
		}
		else {
			lines.push_back(line);
			line = word;
		}
	}
	if (!line.empty())
		lines.push_back(line);
	if (lines.empty())
		lines.emplace_back("");
	return lines;
}

/*
 * Rewrite argv into the form cxxopts accepts:
 *   - any number of leading dashes is reduced to one for single-character
 *     names and two otherwise i.e. '-ref' -> '--ref' and '--N' -> '-N'
 *   - underscores in option names become dashes i.e. '--paired_out' -> '--paired-out'
 *   - a value given as a separate token is attached with '=' so that cxxopts
 *     never has to decide whether the next token is a value. That keeps values
 *     looking like options intact e.g. '--mismatch -3' -> '--mismatch=-3'
 * Tokens that do not name a known option are passed through so that cxxopts
 * reports them, either as an unrecognised option or as an unmatched argument.
 */
std::vector<std::string> Runopts::normalize_argv(int argc, char** argv) const
{
	std::vector<std::string> out;
	out.reserve(static_cast<std::size_t>(argc));
	if (argc > 0)
		out.emplace_back(argv[0]);

	// Split a token into an option descriptor and, if present, an '=' value.
	auto as_option = [this](const std::string& tok, std::string& value, bool& has_value) -> const OptDesc* {
		has_value = false;
		value.clear();
		if (tok.size() < 2 || tok[0] != '-')
			return nullptr;
		auto pos = tok.find_first_not_of('-');
		if (pos == std::string::npos)
			return nullptr; // '-', '--', '---' ...
		std::string body = tok.substr(pos);
		auto eq = body.find('=');
		if (eq != std::string::npos) {
			has_value = true;
			value = body.substr(eq + 1);
			body = body.substr(0, eq);
		}
		return find_opt(body);
	};

	// Is the token an option rather than a value? True for a known option name,
	// and for a dashed multi-character name i.e. a misspelled option. A dashed
	// single character that is not an option is a negative number e.g. '-3'.
	auto looks_like_option = [this](const std::string& tok) -> bool {
		if (tok.size() < 2 || tok[0] != '-')
			return false;
		auto pos = tok.find_first_not_of('-');
		if (pos == std::string::npos)
			return true; // '--', '---' ...
		std::string body = tok.substr(pos);
		auto eq = body.find('=');
		if (eq != std::string::npos)
			body = body.substr(0, eq);
		return find_opt(body) != nullptr || body.size() > 1;
	};

	for (int i = 1; i < argc; ++i) {
		std::string tok(argv[i]);
		if (tok == "--") { // end of options - pass the remainder through verbatim
			for (; i < argc; ++i)
				out.emplace_back(argv[i]);
			break;
		}

		std::string value;
		bool has_value = false;
		const OptDesc* desc = as_option(tok, value, has_value);

		if (desc == nullptr) {
			auto pos = tok.find_first_not_of('-');
			if (tok.size() > 1 && tok[0] == '-' && pos != std::string::npos) {
				// an unrecognised option name. Normalize only the dash prefix so
				// that cxxopts names it in full in the resulting error.
				std::string body = tok.substr(pos);
				out.push_back((body.size() > 1 ? "--" : "-") + body);
			}
			else {
				out.push_back(tok); // a value or a stray argument
			}
			continue;
		}

		const std::string pfx = desc->name.size() > 1 ? "--" : "-";
		if (!has_value && desc->arg != OptArg::NONE && i + 1 < argc
			&& !looks_like_option(argv[i + 1]))
		{
			// the following token is this option's value
			value = argv[i + 1];
			has_value = true;
			++i;
		}
		out.push_back(has_value ? pfx + desc->name + "=" + value : pfx + desc->name);
	}
	return out;
}

/*
 * 'sortmerna <option> -h' | 'sortmerna -h <option>' - help on a single option.
 * Recognized when a help flag is given together with exactly one other token
 * that names an option. The name may be bare ('index') or dashed ('--index').
 */
bool Runopts::print_single_option_help(int argc, char** argv)
{
	if (argc != 3)
		return false;

	const OptDesc* subject = nullptr;
	bool has_help_flag = false;

	for (int i = 1; i < argc; ++i) {
		std::string tok(argv[i]);
		if (tok.empty())
			return false;
		bool is_dashed = tok[0] == '-';
		auto pos = tok.find_first_not_of('-');
		std::string body = (is_dashed && pos != std::string::npos) ? tok.substr(pos) : tok;
		auto canon = canonicalize_opt_name(body);
		if (is_dashed && (canon == OPT_H || canon == OPT_HELP)) {
			has_help_flag = true;
			continue;
		}
		const OptDesc* desc = find_opt(body);
		if (desc == nullptr)
			return false;
		subject = desc;
	}

	if (!has_help_flag || subject == nullptr)
		return false;

	print_option_help(*subject);
	return true;
}

/*
 * '-h' - the compact option list.
 */
void Runopts::print_help()
{
	const std::size_t indent = 4;
	const std::size_t name_w = 20;
	const std::size_t type_w = 14;
	const std::size_t descr_w = 60;
	const std::size_t descr_col = indent + name_w + type_w;

	std::stringstream ss;
	ss << help_header << std::endl;

	std::string cur_group;
	for (const auto& opt : options) {
		std::string group = opt.required ? "[REQUIRED]" : "[" + opt_category_name_map.at(opt.cat) + "]";
		if (group != cur_group) {
			if (!cur_group.empty())
				ss << std::endl;
			ss << std::string(indent, ' ') << group << std::endl;
			cur_group = group;
		}

		std::string name = (opt.name.size() > 1 ? "--" : "-") + opt.name;
		std::string descr = opt.help;
		if (!opt.def.empty())
			descr += " [" + opt.def + "]";

		auto lines = wrap_text(descr, descr_w);
		ss << std::string(indent, ' ')
			<< std::left << std::setw(static_cast<int>(name_w)) << name
			<< std::setw(static_cast<int>(type_w)) << opt.type
			<< lines[0] << std::endl;
		for (std::size_t i = 1; i < lines.size(); ++i)
			ss << std::string(descr_col, ' ') << lines[i] << std::endl;
	}
	std::cout << ss.str();
}

/*
 * '--help' - the full reference, including the extended descriptions.
 */
void Runopts::print_help_full()
{
	const std::size_t indent = 4;
	const std::size_t name_w = 20;
	const std::size_t type_w = 14;
	const std::size_t req_w = 10;
	const std::size_t descr_w = 56;
	const std::size_t descr_col = indent + name_w + type_w + req_w;
	const std::size_t body_indent = indent + 3;

	std::stringstream ss;
	ss << help_header << std::endl;

	std::string cur_group;
	for (const auto& opt : options) {
		std::string group = opt.required ? "[REQUIRED]" : "[" + opt_category_name_map.at(opt.cat) + "]";
		if (group != cur_group) {
			ss << std::endl << std::string(indent, ' ') << group << std::endl << std::endl;
			cur_group = group;
		}

		std::string name = (opt.name.size() > 1 ? "--" : "-") + opt.name;
		auto lines = wrap_text(opt.help, descr_w);
		ss << std::string(indent, ' ')
			<< std::left << std::setw(static_cast<int>(name_w)) << name
			<< std::setw(static_cast<int>(type_w)) << opt.type
			<< std::setw(static_cast<int>(req_w)) << (opt.required ? "Required" : "Optional")
			<< lines[0] << std::endl;
		for (std::size_t i = 1; i < lines.size(); ++i)
			ss << std::string(descr_col, ' ') << lines[i] << std::endl;

		if (!opt.def.empty())
			ss << std::string(descr_col, ' ') << "default: " << opt.def << std::endl;

		if (!opt.help_long.empty()) {
			ss << std::endl;
			std::istringstream body(opt.help_long);
			std::string line;
			while (std::getline(body, line))
				ss << std::string(body_indent, ' ') << line << std::endl;
		}
		ss << std::endl;
	}
	std::cout << ss.str();
}

/*
 * 'sortmerna <option> -h' - everything known about one option.
 */
void Runopts::print_option_help(const OptDesc& opt)
{
	std::string name = (opt.name.size() > 1 ? "--" : "-") + opt.name;
	// the underscore spelling is accepted too, but is not the canonical one
	std::string alias(opt.name);
	std::replace(alias.begin(), alias.end(), '-', '_');

	std::stringstream ss;
	ss << std::endl;
	ss << "  Option:    " << name << std::endl;
	if (alias != opt.name)
		ss << "  Also:      " << (opt.name.size() > 1 ? "--" : "-") << alias
		   << "  (and the single-dash forms '-" << opt.name << "', '-" << alias << "')" << std::endl;
	else if (opt.name.size() > 1)
		ss << "  Also:      -" << opt.name << "  (one or two leading dashes)" << std::endl;
	ss << "  Type:      " << opt.type;
	switch (opt.arg) {
	case OptArg::NONE:     ss << "  (flag, takes no value)"; break;
	case OptArg::REQUIRED: ss << "  (value required)"; break;
	case OptArg::OPTIONAL: ss << "  (value optional)"; break;
	}
	ss << std::endl;
	ss << "  Category:  " << (opt.required ? "REQUIRED" : opt_category_name_map.at(opt.cat)) << std::endl;
	if (!opt.def.empty())
		ss << "  Default:   " << opt.def << std::endl;
	ss << std::endl;

	for (const auto& line : wrap_text(opt.help, 76))
		ss << "  " << line << std::endl;

	if (!opt.help_long.empty()) {
		ss << std::endl;
		std::istringstream body(opt.help_long);
		std::string line;
		while (std::getline(body, line))
			ss << "  " << line << std::endl;
	}
	ss << std::endl;
	std::cout << ss.str();
}

  /*! @fn welcome()
  *  @brief outputs copyright, disclaimer and contact information
  *  @param none
  #  @return none
  */
void about()
{
	std::stringstream ss;

	ss << std::endl
		<< "  Program:      SortMeRNA version " << SORTMERNA_MAJOR << "." << SORTMERNA_MINOR << "." << SORTMERNA_PATCH << std::endl
		<< "  Description:  Filters ribosomal RNA (rRNA) reads from metagenomic and" << std::endl
		<< "                metatranscriptomic sequencing datasets." << std::endl
		<< "  Copyright:    2016-2026 Clarity Genomics Inc" << std::endl
		<< "                2014-2016 Knight Lab, UCSD" << std::endl
		<< "                2012-2014 Bonsai Bioinformatics Research Group, University Lille 1 / Inria Nord-Europe" << std::endl
		<< "  Disclaimer:   SortMeRNA comes with ABSOLUTELY NO WARRANTY; without even the" << std::endl
		<< "                implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl
		<< "                See the GNU Lesser General Public License for more details." << std::endl;

	std::cout << ss.str();
}
} // namespace sortmerna
