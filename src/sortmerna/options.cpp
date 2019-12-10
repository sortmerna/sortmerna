/* 
 * FILE: options.cpp
 * Created: Jun 07, 2018 Thu
 * @copyright 2016-19 Clarity Genomics BVBA
 */

 // TODO: BUG: if SMR headers moved down after 3rd party, 'timeval' struct gets 'redefined' - compiler error. That's a header mess bug.

 // 3rd party
#include "zlib.h"

#include "version.h"
#include "build_version.h"
#include "options.hpp"
#include "common.hpp"
#include "gzip.hpp"
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


#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

/*! @brief Measure time using this variable. */
//timeval t;

// forward
void about();
//void help();
std::string get_user_home(); // util.cpp
//unsigned int list_dir(std::string dpath); // replaced with filesystem
//bool dirExists(std::string dpath); // replaced with filesystem
std::string trim_leading_dashes(std::string const& name); // util.cpp
std::string get_basename(const std::string &file); // util.cpp
std::streampos filesize(const std::string &file); // util.cpp
std::string string_hash(const std::string &val); // util.cpp

Runopts::Runopts(int argc, char**argv, bool dryrun)
{
	process(argc, argv, dryrun);
	if (skiplengths.empty())
	{
		for (int i = 0; i < indexfiles.size(); ++i)
		{
			skiplengths.push_back({ 0,0,0 });
		}
	}
}

/* 
 * validates file (be it plain or zipped) and sets corresponding options
 * For paired reads files use --reads file1 --reads file2
 */
void Runopts::opt_reads(const std::string &file)
{
	std::stringstream ss;

	auto numread = mopt.count("reads");
	auto readcnt = readfiles.size();

	std::cout << STAMP << "Processing reads file [" << readcnt + 1 << "] out of total [" << numread << "] files" << std::endl;

	if (file.size() == 0)
	{
		ERR(help_reads);
		exit(EXIT_FAILURE);
	}

	// check file exists and can be read
	auto fpath = std::filesystem::path(file);
	auto fpath_a = std::filesystem::path(); // absolute path

	if (std::filesystem::exists(fpath))
	{
		fpath_a = fpath;
	}
	else if (fpath.is_relative())
	{
		std::cout << "File  [" << file << "] appears to be a relative path. Trying to locate in current directory and in the working directory ...";
		fpath_a = std::filesystem::current_path() / file; // std::filesystem::path(file) makes no difference
		if (std::filesystem::exists(fpath_a))
		{
			std::cout << "File [" << fpath_a << "] exists" << std::endl;
		}
		else
		{
			std::cout << "File  [" << file << "] has not been found in the current directory. Trying the working directory ...";
			fpath_a = std::filesystem::path(workdir) / file;
			if (!std::filesystem::exists(fpath_a))
			{
				std::cout << STAMP
					<< "Could not locate File [" << file << "] neither at current path [" << std::filesystem::current_path()
					<< "] nor in Workdir [" << workdir << "]" << std::endl;
				ERR(help_reads);
				exit(EXIT_FAILURE);
			}
		}
	}
	else
	{
		ss << STAMP <<
			"The file [" << file << "] is not an existing/valid absolute or relative path\n" << help_reads;
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	std::ifstream ifs(fpath_a, std::ios_base::in | std::ios_base::binary);
	if (!ifs.is_open())
	{
		ss << STAMP << "Failed to open file [" << fpath_a << "]";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	bool gzipped = "gz" == file.substr(file.rfind('.') + 1); // file ends with 'gz'
	std::string line;
	Gzip gzip(gzipped);
	int stat = gzip.getline(ifs, line);
	if (RL_OK == stat)
	{
		is_gz = gzipped;
	}
	else if (!gzipped)
	{
		// try reading as gzipped even though file has no 'gz' extension
		gzip.~Gzip(); // destroy
		Gzip gzip(true);
		stat = gzip.getline(ifs, line);
		if (RL_OK == stat) is_gz = true;
	}

	if (RL_OK == stat && line.size() > 0)
	{
		have_reads = true;
		readfiles.push_back(fpath_a.generic_string());
	}
	else
	{
		ss << "Could not read from file " << file << " [" << strerror(errno) << "]";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	if (ifs.is_open())
		ifs.close();
} // ~Runopts::opt_reads

void Runopts::opt_ref(const std::string &file)
{
	std::stringstream ss;

	auto numref = mopt.count("ref");
	auto refcnt = indexfiles.size();

	std::cout << STAMP << "Processing reference [" << refcnt + 1 << "] out of total [" << numref << "] references" << std::endl;

	if (file.size() == 0)
	{
		ERR(help_ref);
		exit(EXIT_FAILURE);
	}

	// check file exists and can be read
	auto fpath = std::filesystem::path(file);
	auto fpath_a = std::filesystem::path(); // absolute path

	if (std::filesystem::exists(fpath))
	{
		fpath_a = fpath;
	}
	else if (fpath.is_relative())
	{
		std::cout << "File  [" << file << "] appears to be a relative path. Trying to locate in current directory and in the working directory ...";
		fpath_a = std::filesystem::current_path() / file; // std::filesystem::path(file) makes no difference
		if (std::filesystem::exists(fpath_a))
		{
			std::cout << "File [" << fpath_a << "] exists" << std::endl;
		}
		else
		{
			std::cout << "File  [" << file << "] has not been found in the current directory. Trying the working directory ...";
			fpath_a = std::filesystem::path(workdir) / file;
			if (!std::filesystem::exists(fpath_a))
			{
				std::cout << STAMP 
					<< "Could not locate File [" << file << "] neither at current path [" << std::filesystem::current_path() 
					<< "] nor in Workdir [" << workdir << "]" << std::endl;
				ERR(help_ref);
				exit(EXIT_FAILURE);
			}
		}
	}
	else
	{
		ss << STAMP
			<< "The file [" << file << "] is not an existing/valid absolute or relative path\n"
			<< help_ref;
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// check files are readable
	std::ifstream ifstr(fpath_a);
	if (!ifstr.is_open() || !ifstr.good()) {
		ss << STAMP << "Cannot read file [" << file << "]";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}
	else {
		ss << STAMP << "File [" << std::filesystem::absolute(file) << "] exists and is readable\n";
		std::cout << ss.str();
	}

	if (ifstr.is_open())
		ifstr.close();

	// check index file names are distinct
	for (int i = 0; i < (int)indexfiles.size(); i++)
	{
		if ((indexfiles[i].first).compare(file) == 0)
		{
			ss.str("");
			ss << STAMP << "Reference file (" << file << ") has been entered more than once. Ignoring redundant enties";
			WARN(ss.str());
		}
	}

	// derive index file prefix from the reference file name
	// if we are here the Workdir is OK
	//std::string basename = get_basename(file);
	//std::string idx_file_pfx = workdir + "/" + IDX_DIR + "/" + string_hash(basename);
	auto fpath_base = std::filesystem::path(file).filename();
	auto idx_file_pfx = std::filesystem::path(workdir) / IDX_DIR / string_hash(fpath_base.generic_string());

	indexfiles.push_back(std::pair<std::string, std::string>(fpath_a.generic_string(), idx_file_pfx.generic_string()));
} // ~Runopts::opt_ref

void Runopts::opt_aligned(const std::string &file)
{
	if (file.size() == 0)
	{
		std::cout << STAMP << "File name was not provided with option '--aligned [FILE]'. Using default name 'aligned'" << std::endl;
	}
} // ~Runopts::opt_aligned

/**
 * depends on '--fastxout'
 */
void Runopts::opt_other(const std::string &file)
{
	auto cnt = mopt.count("fastxout");
	if (cnt < 0)
	{
		ERR("Option '--other' can only be used together with '--fastx' option.");
		exit(EXIT_FAILURE);
	}

	if (file.size() == 0)
	{
		std::cout << STAMP << "'--other' option was specified. Using file name: [" << other_out_pfx << "]" << std::endl;
	}
	is_other = true;
} // ~Runopts::opt_other

void Runopts::opt_log(const std::string &val)
{
	if (is_log)
		WARN("'--" << OPT_LOG << "' is deprecated. True by default.");
} // ~Runopts::optLog

void Runopts::opt_de_novo_otu(const std::string &val)
{
	is_de_novo_otu = true;
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

void Runopts::opt_paired_in(const std::string &val)
{
	is_paired_in = true;
} // ~Runopts::optPairedIn

void Runopts::opt_paired_out(const std::string &val)
{
	is_paired_out = true;
} // ~Runopts::optPairedOut

void Runopts::opt_match(const std::string &val)
{
	std::stringstream ss;
	if (val.size() == 0)
	{
		ss << STAMP << "'" << OPT_MATCH << "' " << "requires a positive integer as input e.g. 2";
		ERR(ss.str());
		exit(EXIT_FAILURE);
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
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}
} // ~Runopts::optMatch

void Runopts::opt_mismatch(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR("--mismatch [INT] requires a negative integer input (ex. --mismatch -2)");
		exit(EXIT_FAILURE);
	}

	// set mismatch
	if (!mismatch_set)
	{
		mismatch = atoi(val.data());
		if (mismatch > 0)
		{
			ERR("--mismatch [INT] takes a negative integer (ex. --mismatch -2)");
			exit(EXIT_FAILURE);
		}
		mismatch_set = true;
	}
	else
	{
		ERR("--mismatch [INT] has been set twice, please verify your choice");
		print_help();
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_mismatch

void Runopts::opt_gap_open(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR("--gap_open [INT] requires a positive integer as input (ex. --gap_open 5)");
		exit(EXIT_FAILURE);
	}

	// set gap open
	if (!gap_open_set)
	{
		gap_open = atoi(val.data());
		if (gap_open < 0)
		{
			ERR("--gap_open [INT] requires a positive integer as input (ex. --gap_open 5)");
			exit(EXIT_FAILURE);
		}
		gap_open_set = true;
	}
	else
	{
		ERR("--gap_open [INT] has been set twice, please verify your choice");
		print_help();
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_gap_open

void Runopts::opt_gap_ext(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR("--gap_ext [INT] requires a positive integer as input (ex. --gap_ext 2)");
		exit(EXIT_FAILURE);
	}
	// set gap extend
	if (!gap_ext_set)
	{
		gap_extension = atoi(val.data());
		if (gap_extension < 0)
		{
			ERR("--gap_ext [INT] requires a positive integer as input (ex. --gap_ext 2)");
			exit(EXIT_FAILURE);
		}
		gap_ext_set = true;
	}
	else
	{
		ERR("--gap_ext [INT] has been set twice, please verify your choice");
		print_help();
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_gap_ext

void Runopts::opt_num_seeds(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR("--num_seeds [INT] requires a positive integer as input (ex. --num_seeds 6)");
		exit(EXIT_FAILURE);
	}
	// set number of seeds
	if (seed_hits < 0)
	{
		char* end = 0;
		seed_hits = (int)strtol(val.data(), &end, 10); // convert to integer
		if (seed_hits <= 0)
		{
			ERR("--num_seeds [INT] requires a positive integer (>0) as input (ex. --num_seeds 6)");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		ERR("--num_seeds [INT] has been set twice, please verify your choice");
		print_help();
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_num_seeds

/* --fastx */
void Runopts::opt_fastx(const std::string &val)
{
	if (is_fastx)
	{
		ERR("--fastx has already been set once.");
		exit(EXIT_FAILURE);
	}
	else
	{
		is_fastx = true;
	}
} // ~Runopts::opt_fastx

void Runopts::opt_sam(const std::string &val)
{
	std::stringstream ss;
	if (is_sam)
	{
		ss << STAMP << "'" << OPT_SAM << "' has already been set. Ignoring second flag.";
		WARN(ss.str());
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
			ss.str("");
			ss << ": `" << opt << "` is not supported in --blast [STRING].";
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
	}
	// more than 1 field with blast human-readable format given
	if (blast_human_readable && (blastops.size() > 1))
	{
		ss.str("");
		ss << ": for human-readable format, --blast [STRING] can only contain a single field '0'.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}
	// both human-readable and tabular format options have been chosen
	if (blast_human_readable && blastFormat == BlastFormat::TABULAR)
	{
		ss.str("");
		ss << ": --blast [STRING] can only have one of the options '0' (human-readable) or '1' (tabular).";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	is_blast = true;
} // ~Runopts::opt_blast

void Runopts::opt_min_lis(const std::string &val)
{
	std::stringstream ss;
	if (val.size() == 0)
	{
		ss << STAMP << "'" << OPT_MIN_LIS
			<< "' [INT] requires a positive integer as input e.g. 2. (if 0 - all high scoring reference sequences are searched)";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// min_lis_gv has already been set
	if (is_min_lis)
	{
		ss.str("");
		ss << STAMP << "'" << OPT_MIN_LIS << "' [INT] has been set twice, please verify your choice.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}
	else
	{
		if ((sscanf(val.data(), "%d", &min_lis) != 1) || (min_lis < 0))
		{
			ss.str("");
			ss << STAMP << "'" << OPT_MIN_LIS
				<< "' [INT] requires a positive integer as input e.g. 2. If 0, all high scoring reference sequences are searched)";
			exit(EXIT_FAILURE);
		}
		is_min_lis = true;
	}
} // ~Runopts::opt_min_lis

void Runopts::opt_best(const std::string &val)
{
	std::stringstream ss;
	if (val.size() == 0)
	{
		ss << STAMP << "'" << OPT_BEST << "'" << " [INT] requires a positive integer e.g. 2).";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	if (is_best)
	{
		ss.str("");
		ss << STAMP << "'" << OPT_BEST << "'" << " [INT] has been set twice, please verify your choice.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}
	else
	{
		if ((sscanf(val.data(), "%d", &num_best_hits) != 1))
		{
			ss.str("");
			ss << STAMP << "Could not read '" << OPT_BEST << "' [INT] as integer";
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
		is_best = true;
	}
} // ~Runopts::opt_best

void Runopts::opt_num_alignments(const std::string &val)
{
	std::stringstream ss;
	if (val.size() == 0)
	{
		ss << STAMP << "'" << OPT_NUM_ALIGNMENTS
			<< "' [INT] requires a posistive integer as input e.g. 2. If 0, all alignments are output.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	if (is_num_alignments)
	{
		ss.str("");
		ss << STAMP << "'" << OPT_NUM_ALIGNMENTS 
			<< "' [INT] has been set twice, please verify your parameters.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// set number of alignments to output reaching the E-value
	num_alignments = atoi(val.data());
	if (num_alignments < 0)
	{
		ss.str("");
		ss << STAMP << "'" << OPT_NUM_ALIGNMENTS
			<< "' [INT] requires a posistive integer as input e.g. 2. If 0, all alignments are output.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}
	is_num_alignments = true;
} // ~Runopts::opt_num_alignments

void Runopts::opt_edges(const std::string &val)
{
	// --edges is already set
	if (edges_set)
	{
		ERR(" : --edges [INT]%% has already been set once.");
		exit(EXIT_FAILURE);
	}

	char *end = 0;
	// find if % sign exists
	if (val.find_first_of("%") != std::string::npos)
		is_as_percent = true;

	// convert to integer
	edges = std::stoi(val); //edges = (int)strtol(val.data(), &end, 10);

	if (edges < 1 || edges > 10)
	{
		ERR(" : --edges [INT]%% requires a positive integer between 0-10 as input (ex. --edges 4).");
		exit(EXIT_FAILURE);
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
		ERR(" : BOOL --SQ has been set twice, please verify your choice.");
		exit(EXIT_FAILURE);
	}

	is_SQ = true;
} // ~Runopts::optSQ

void Runopts::opt_passes(const std::string &val)
{
	std::stringstream ss;
	if (passes_set)
	{
		ss << STAMP << "'" << OPT_PASSES << "' [INT,INT,INT] has been set twice, please verify your choice.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// set passes
	for (auto pos = val.find(",", 0), pos_from = pos-pos, count = pos-pos; pos != std::string::npos; pos = val.find(",", pos_from))
	{
		auto tok = val.substr(pos_from, pos);
		pos_from += pos +1;
		if (++count > 3) 
		{
			ss.str("");
			ss << STAMP << "Exactly 3 integers has to be provided with '" << OPT_PASSES << "' [INT,INT,INT]";
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
		auto skiplen = std::stoi(tok);
		if (skiplen > 0)
			skiplengths.emplace_back(skiplen);
		else
		{
			ss.str("");
			ss << STAMP << "All three integers in '" << OPT_PASSES
				<< "' [INT,INT,INT] must contain positive integers where 0 < INT < (shortest read length).";
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}
	}
	passes_set = true;
} // ~Runopts::opt_passes

void Runopts::opt_id(const std::string &val)
{
	// % id
	if (align_id < 0)
	{
		if ((sscanf(val.data(), "%lf", &align_id) != 1) ||
			(align_id < 0) || (align_id > 1))
		{
			ERR(" : --id [DOUBLE] must be a positive float with value 0<=id<=1.");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		ERR(" : --id [DOUBLE] has been set twice, please verify your command parameters.");
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_id

void Runopts::opt_coverage(const std::string &val)
{
	// % query coverage
	if (align_cov < 0)
	{
		if ((sscanf(val.data(), "%lf", &align_cov) != 1) ||
			(align_cov < 0) || (align_cov > 1))
		{
			ERR(" : --coverage [DOUBLE] must be a positive float with value 0<=id<=1.");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		ERR(" : --coverage [DOUBLE] has been set twice, please verify your command parameters.");
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_coverage

void Runopts::opt_version(const std::string &val)
{
	std::cout << std::endl
		<< "SortMeRNA version " << SORTMERNA_MAJOR << "." << SORTMERNA_MINOR << "." << SORTMERNA_PATCH << std::endl
		<< "Build Date: " << sortmerna_build_compile_date << std::endl
		<< sortmerna_build_git_sha << std::endl
		<< sortmerna_build_git_date << std::endl;
	exit(EXIT_SUCCESS);
} // ~Runopts::opt_version

void Runopts::opt_unknown(char **argv, int &narg, char * opt)
{
	std::stringstream ss;
	ss << " : option --" << opt << " not recognized";
	ERR(ss.str());
	print_help();
	exit(EXIT_FAILURE);
} // ~Runopts::opt_unknown

void Runopts::opt_e(const std::string &val)
{
	std::stringstream ss;

	// E-value
	if (val.size() == 0)
	{
		ERR(": -e [DOUBLE] requires a positive double as input (ex. --e 1e-5)");
		exit(EXIT_FAILURE);
	}

	if (evalue < 0)
	{
		sscanf(val.data(), "%lf", &evalue);
		if (evalue < 0)
		{
			ERR(" : -e [DOUBLE] requires a positive double as input (ex. --e 1e-5)");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		ERR(" : -e [DOUBLE] has been set twice, please verify your command parameters.");
		exit(EXIT_FAILURE);
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
		ERR(" : BOOL -F has been set more than once, please check your command parameters.");
		exit(EXIT_FAILURE);
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
		ERR(" : BOOL '-R' has been set more than once, please check your command parameters.");
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_R

void Runopts::opt_h(const std::string &val)
{
	about();
	print_help();
	exit(0);
} // ~Runopts::opt_h

void Runopts::opt_v(const std::string &val)
{
	is_verbose = true;
} // ~Runopts::opt_v

void Runopts::opt_N(const std::string &val)
{
	if (!match_ambiguous_N)
	{
		match_ambiguous_N = true;
		score_N = std::stoi(val);
	}
	else
	{
		ERR(": BOOL -N has been set more than once, please check your command parameters.");
		exit(EXIT_FAILURE);
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
		ERR(msg);
		exit(EXIT_FAILURE);
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


void Runopts::opt_thpp(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR(": --thpp [INT:INT] requires 2 integers for number of "
			<< "Read:Processor threads (ex. --thpp 1:1)");
		exit(EXIT_FAILURE);
	}

	std::istringstream strm(val);
	std::string tok;
	for (int i = 0; std::getline(strm, tok, ':'); ++i)
	{
		switch (i)
		{
		case 0: num_read_thread_pp = std::stoi(tok); break;
		case 1: num_proc_thread_pp = std::stoi(tok); break;
		}
	}
} // ~Runopts::opt_threads_pp


void Runopts::opt_threp(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR(": --threp [INT:INT] requires 2 integers for number of "
			<< "Read:Processor threads (ex. --threp 1:1)");
		exit(EXIT_FAILURE);
	}

	std::istringstream strm(val);
	std::string tok;
	for (int i = 0; std::getline(strm, tok, ':'); ++i)
	{
		switch (i)
		{
		case 0: num_read_thread_rep = std::stoi(tok); break;
		case 1: num_proc_thread_rep = std::stoi(tok); break;
		}
	}
} // ~Runopts::opt_threads

// KeyValDatabase
void Runopts::opt_d(const std::string &val)
{
	if (val.size() != 0)
		kvdbPath = val;
} // ~Runopts::opt_d


void Runopts::opt_dbg_put_db(const std::string &val)
{
	is_dbg_put_kvdb = true;
}

void Runopts::opt_default(const std::string &opt)
{
	std::stringstream ss;
	ss << STAMP << "Option: '" << opt << "' is not recognized";
	ERR(ss.str());
	print_help();
	exit(EXIT_FAILURE);
} // ~Runopts::opt_default

  /* Processing task */
void Runopts::opt_task(const std::string &val)
{
	int task_num = std::stoi(val);

	if (task_num > 4)
	{
		std::stringstream ss;
		ss << "Option '" << OPT_TASK << "' can only take values in range [0..4] Provided value is [" << task_num << "'";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	switch (task_num)
	{
	case 0: alirep = align; break;
	case 1: alirep = postproc; break;
	case 2: alirep = report; break;
	case 3: alirep = alipost; break;
	case 4: alirep = all; break;
	}
} // ~Runopts::optReport

// interactive session '--cmd'
void Runopts::opt_cmd(const std::string &val)
{
	is_cmd = true;
} // ~Runopts::optInteractive

void Runopts::opt_workdir(const std::string &path)
{
	if (path.size() == 0)
	{
		auto wdpath = std::filesystem::path(get_user_home()) / WORKDIR_DEF_SFX;
		workdir = wdpath.string();
		std::cout << STAMP << "'" << OPT_WORKDIR 
			<< "' option was not provided. Using USERDIR to set the working directory: [" << workdir << "]" << std::endl;
	}
	else {
		workdir = path;
		std::cout << STAMP << "Using WORKDIR [" << std::filesystem::absolute(workdir) << " as specified" << std::endl;
	}
}

// indexing options
void Runopts::opt_tmpdir(const std::string &val)
{
	std::cout << STAMP << "TODO: deprecated indexing option: " << help_tmpdir << " To be removed." << std::endl;
}

void Runopts::opt_interval(const std::string &val)
{
	std::stringstream ss;
	auto count = mopt.count("interval");
	if (count > 1)
	{
		ss << " Option 'interval' entered [" << count << "] times. Only the last value will be used. " << help_interval;
		WARN(ss.str());
	}

	if (val.size() == 0)
	{
		WARN("Option 'interval' given without value. Using default: " + interval);
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
	std::stringstream ss;
	auto count = mopt.count("m");
	if (count > 1)
	{
		ss << " Option 'm' entered [" << count << "] times. Only the last value will be used. " << help_m;
		WARN(ss.str());
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
	std::stringstream ss;
	auto count = mopt.count(OPT_L);
	if (count > 1)
	{
		ss << " Option '" << OPT_L << "' entered [" << count << "] times. Only the last value will be used." << std::endl << "\tHelp: " << help_L;
		WARN(ss.str());
	}

	if (val.size() == 0)
	{
		ss.str("");
		ss << "Option 'L' given without value. Using default: " << seed_win_len;
		WARN(ss.str());
	}
	else
	{
		int lnwin_t = std::stoi(val);

		if (lnwin_t <= 0 || lnwin_t % 2 == 1 || lnwin_t < 8 || lnwin_t > 26)
		{
			ss.str("");
			ss << STAMP 
				<< "Option '" << OPT_L <<"' takes a Positive Even integer between 8 and 26 inclusive e.g. 10, 12, 14, .. , 20. Provided value: " 
				<< lnwin_t << " Default will be used: " << seed_win_len;
			WARN(ss.str());
		}
		else
		{
			seed_win_len = lnwin_t;
		}
	}
} // ~Runopts::opt_L

void Runopts::opt_max_pos(const std::string &val)
{
	std::stringstream ss;
	auto count = mopt.count(OPT_MAX_POS);
	if (count > 1)
	{
		ss << " Option '" << OPT_MAX_POS << "' entered [" << count << "] times. Only the last value will be used" << std::endl 
			<< "\tHelp: " << help_max_pos;
		WARN(ss.str());
	}

	if (val.size() == 0)
	{
		ss.str("");
		ss << "Options 'max_pos' takes a positive integer e.g. 250. Using default: " << max_pos;
		WARN(ss.str());
	}
	else
	{
		max_pos = std::stoi(val);
	}
}

void Runopts::test_kvdb_path()
{
	if (kvdbPath.size() == 0)
	{
		kvdbPath = workdir + "/" + KVDB_DIR;
	}

	std::cout << STAMP << "Key-value DB location (" << std::filesystem::absolute(kvdbPath) << ")" << std::endl;

	if (std::filesystem::exists(kvdbPath))
	{
		// dir exists and is empty
		//auto count = list_dir(kvdbPath);
		if (std::filesystem::is_empty(kvdbPath))
		{
			// dir exists and empty -> use
			if (ALIGN_REPORT::postproc == alirep || ALIGN_REPORT::report == alirep)
			{
				std::cout << STAMP << "KVDB (" << std::filesystem::absolute(kvdbPath) << " is empty. Will run alignment" << std::endl;
			}
		}
		else // not empty
		{
			// TODO: Store some metadata in DB to verify the alignment.
			// kvdb.verify()
			if (ALIGN_REPORT::align == alirep || ALIGN_REPORT::all == alirep || ALIGN_REPORT::alipost == alirep)
			{
				// if (kvdb.verify()) // TODO
				// output the listing
				std::stringstream ss;
				ss << STAMP << "Path '" << std::filesystem::absolute(kvdbPath) << "' exists with the following content:" << std::endl;

				for (auto& subpath : std::filesystem::directory_iterator(kvdbPath))
					ss << subpath.path().filename() << std::endl;

				ss << "\tPlease, ensure the directory [" << std::filesystem::absolute(kvdbPath) << "] is Empty prior running 'sortmerna'" << std::endl;
				WARN(ss.str());
				exit(EXIT_FAILURE);
			}
		}
	}
	else
	{
		// dir does not exist -> try creating
		std::cout << STAMP << "Database (" << std::filesystem::absolute(kvdbPath) << ") will be created" << std::endl;
	}
} // ~test_kvdb_path

/** 
 * main method of this class. 
 * Parses the command line options, validates, and sets the class member variables 
 */
void Runopts::process(int argc, char**argv, bool dryrun)
{
#if defined(__APPLE__)
	int sz[2] = { CTL_HW, HW_MEMSIZE };
	u_int namelen = sizeof(sz) / sizeof(sz[0]);
	uint64_t size;
	size_t len = sizeof(size);
	if (sysctl(sz, namelen, &size, &len, NULL, 0) < 0)
	{
		fprintf(stderr, "\n  %sERROR%s: sysctl (main.cpp)\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
#endif

#if defined(_WIN32)
	_setmode(_fileno(stderr), _O_BINARY);
#endif

	if (argc == 1)
	{
		about();
		ERR("Missing required command options");
		exit(EXIT_FAILURE);
	}

	std::cout << "\n" << STAMP << "=== Options processing starts ... ===\n\n";

	// store the command line
	for (int i = 0; i < argc; i++)
	{
		cmdline.append(argv[i]);
		cmdline.append(" ");
	}

	std::string flag;

	// parse cmd options and store into 'mopt' multimap
	for (auto i = argc - argc, flag_count = 0; i != argc; ++i)
	{
		// if arg starts with dash it is flag
		bool is_flag = '-' == **(argv + i); // first character is dash '-'
		if (is_flag && flag_count == 0) {
			std::cout << "Found flag: " << *(argv + i) << std::endl;
			if (i == argc - 1) {
				flag = trim_leading_dashes(*(argv + i));
				mopt.emplace(std::make_pair(flag, "")); // add the last boolean flag
			}
			++flag_count;
			continue;
		}

		if (is_flag && flag_count == 1) {
			std::cout << "Previous flag: " << *(argv + i - 1) << " is Boolean. Setting to True" << std::endl;
			std::cout << "Found flag: " << *(argv + i) << std::endl;
			flag = trim_leading_dashes(*(argv + i - 1));
			mopt.emplace(std::make_pair(flag, "")); // previous boolean flag
			if (i == argc - 1) {
				flag = trim_leading_dashes(*(argv + i));
				mopt.emplace(std::make_pair(flag, "")); // add the last boolean flag
			}
			continue;
		}

		if (!is_flag && flag_count == 0) {
			std::cout << "Found value: " << *(argv + i) << std::endl;
			continue;
		}

		if (!is_flag && flag_count == 1) {
			std::cout << "Found value: " << *(argv + i) << " of previous flag: " << *(argv + i - 1) << std::endl;
			flag = trim_leading_dashes(*(argv + i - 1));
			mopt.emplace(std::make_pair(flag, *(argv + i)));
			--flag_count;
			continue;
		}
	} // ~for parsing cmd options

	// is this only Help options?
	bool is_help_opt = mopt.count(OPT_H) + mopt.count(OPT_VERSION) > 0;

	if (!is_help_opt)
	{
		// check required options were provided
		for (auto opt : options)
		{
			if (std::get<3>(opt))
			{
				if (mopt.count(std::get<0>(opt)) == 0)
				{
					std::cout << "Missing required flag: " << std::get<0>(opt) << std::endl;
				}
			}
		}

		// process WORKDIR first as other options depend on it
		auto wd_it = mopt.find("workdir");
		std::string wdir = "";
		if (wd_it != mopt.end())
		{
			wdir = wd_it->second;
			mopt.erase(wd_it); // remove to prevent repeated processing below
		}
		opt_workdir(wdir);
	}

	// loop through the rest of options
	for (auto opt : mopt)
	{
		std::cout << STAMP << "Processing option: " << opt.first << " with value: " << opt.second << std::endl;

		bool has_opt = false;
		for (auto optx : options)
		{
			if (std::get<0>(optx) == opt.first)
			{
				// call processing function for the given option
				std::invoke(std::get<5>(optx), this, opt.second);
				has_opt = true;
			}
		}

		if (!has_opt) 
		{
			// passed option is not recognized
			opt_default(opt.first);
		}
	}

	std::cout << "\n" << STAMP << "=== Options processing done ===\n\n";

	if (!is_help_opt)
	{
		// validate the options. TODO: should be part of options validation?
		test_kvdb_path();
		validate();
		about(); // if we are here, args are OK, welcome the user
	}
} // ~Runopts::process

/**
 * Validate the options and setup defaults
 */
void Runopts::validate()
{
	std::stringstream ss;

	// No output format has been chosen
	if (!(is_fastx || is_blast || is_sam || is_otu_map || is_de_novo_otu))
	{
		is_blast = true;
		std::cout << STAMP 
			<< "No output format has been chosen (fastx/sam/blast/otu_map). Using default '" 
			<< OPT_BLAST << "'" << std::endl;
	}

	if (is_paired_in && is_paired_out)
	{
		ss << STAMP << "Options '" << OPT_PAIRED_IN << "' and '" << OPT_PAIRED_OUT << "' are mutually exclusive. Please choose one or the other";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// Options --paired_in and --paired_out can only be used with FASTA/Q output
	if (!is_fastx && (is_paired_in || is_paired_out))
	{
		ss.str("");
		ss << STAMP << "Options '" << OPT_PAIRED_IN << "' and '" << OPT_PAIRED_OUT
			<< "' must be accompanied by option '" << OPT_FASTX << "'. Setting to true.";
		is_fastx = true;
	}

	// An OTU map can only be constructed with the single best alignment per read
	if (is_otu_map && is_num_alignments)
	{
		ERR("'--otu_map' cannot be set together with --num_alignments [INT].\n"
			"\tThe option --num_alignments [INT] doesn't keep track of"
			"\tthe best alignment which is required for constructing an OTU map.\n"
			"\tUse --otu_map with --best [INT] instead.");
		exit(EXIT_FAILURE);
	}

	if (is_num_alignments && !(is_blast || is_sam || is_fastx))
	{
		ss.str("");
		ss << STAMP << "'" << OPT_NUM_ALIGNMENTS
			<< "' [INT] has been set but no output format has been chosen (--blast | --sam | --fastx). Using default '" 
			<< OPT_BLAST << "'";
		WARN(ss.str());
		exit(EXIT_FAILURE);
	}

	// If --best output was chosen, check an alignment format has also been chosen
	if (is_best && !(is_blast || is_sam || is_otu_map))
	{
		ss.str("");
		ss << STAMP << "'" << OPT_BEST 
			<< "' [INT] has been set but no output format has been chosen (--blast | --sam | --otu_map). Using default '" 
			<< OPT_BLAST << "'";
		WARN(ss.str());
	}

	// Check gap extend score < gap open score
	if (gap_extension > gap_open)
	{
		std::stringstream ss;
		ss << STAMP << "--gap_ext [INT] must be less than --gap_open [INT].";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// Option --print_all_reads can only be used with Blast-like tabular
	// and SAM formats (not pairwise)
	if (is_print_all_reads && is_blast && blastFormat != BlastFormat::TABULAR)
	{
		std::stringstream ss;
		ss << STAMP << "--print_all_reads [BOOL] can only be used with BLAST-like";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// Only one of these options is allowed (--best outputs one alignment,
	// --num_alignments outputs > 1 alignments)
	if (is_best && is_num_alignments)
	{
		std::stringstream ss;
		ss << STAMP 
			<< "'"<< OPT_BEST << "' [INT] and '" 
			<< OPT_NUM_ALIGNMENTS << "' [INT] cannot be set together.\n'" 
			<< OPT_BEST <<"' searches [INT] highest scoring reference sequences\n"
				"and outputs a single best alignment, whereas '" 
			<< OPT_NUM_ALIGNMENTS << "'\n" << "outputs the first [INT] alignments.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// Option --min_lis [INT] can only accompany --best [INT]
	if (is_min_lis && is_num_alignments)
	{
		std::stringstream ss;
		ss << STAMP << "'" << OPT_MIN_LIS << "' [INT] and '" << OPT_NUM_ALIGNMENTS << "' [INT] cannot be set together.\n";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// Option --mis_lis INT accompanies --best INT, cannot be set alone
	if (is_min_lis && !is_best)
	{
		std::stringstream ss;
		ss << STAMP	<< "--min_lis [INT] must be set together with --best [INT].";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// %id and %coverage can only be used with --otu_map
	if ((align_id > 0 || align_cov > 0) && !is_otu_map)
	{
		std::stringstream ss;
		ss << STAMP
			<< "--id [INT] and --coverage [INT] can only be used together with --otu_map.\n"
			"\tThese two options are used for constructing the OTU map\n"
			"\tby filtering alignments passing the E-value threshold.";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	// if neither strand was selected for search, search both
	if (!is_forward && !is_reverse)
	{
		is_forward = true;
		is_reverse = true;
	}
	// default number of threads is 1
	//if (numcpu_gv < 0) numcpu_gv = 1;
	// default E-value 
	if (evalue < 0.0) evalue = 1;

	// SW alignment parameters
	if (!match_set) 
		match = 2;

	if (!mismatch_set) 
		mismatch = -3;

	if (!gap_open_set) 
		gap_open = 5;

	if (!gap_ext_set) 
		gap_extension = 2;

	if (!match_ambiguous_N) 
		score_N = mismatch;

	// default method for searching alignments
	if (!is_best && !is_num_alignments)
	{
		// TODO: looks arbitrary. Why the alignment contolling options would depend on the output?
		// FASTA/FASTQ output, stop searching for alignments after the first match
		if (is_fastx && !(is_blast || is_sam || is_otu_map || is_log || is_de_novo_otu))
			num_alignments = 1;
		// output single best alignment from best candidate hits
		else
		{
			num_best_hits = 1;
			min_lis = 2;
		}
	}

	// default minimum LIS used for setting the number of
	// alignments to search prior to outputting --best INT
	if (is_best && !is_min_lis) 
		min_lis = 2;

	// default number of seed hits before searching for candidate LIS
	if (seed_hits < 0) 
		seed_hits = 2;

	// default number of nucleotides to add to each edge of an alignment
	// region before extension
	if (edges < 0) 
		edges = 4;

	// activate heuristic for stopping search (of 1-error matches) after
	// finding 0-error match
	if (!full_search_set) 
		is_full_search = false;

	// default %id to keep alignment
	if (align_id < 0)
	{
		// if OTU-map is chosen, set default similarity to 0.97
		if (is_otu_map) align_id = 0.97;
		else align_id = 0;
	}

	// default %query coverage to keep alignment
	if (align_cov < 0)
	{
		// if OTU-map is chosen, set default coverage to 0.97
		if (is_otu_map) align_cov = 0.97;
		else align_cov = 0;
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

void Runopts::print_help()
{
	std::stringstream ss;
	std::string req;
	std::string pfx;

	bool set_req = false;
	bool set_com = false;
	bool set_otu = false;
	bool set_adv = false;
	bool set_idx = false;
	bool set_hlp = false;
	bool set_dev = false;

	std::string HDR_REQ = "[REQUIRED]";
	std::string HDR_COM = "[COMMON]";
	std::string HDR_OTU = "[OTU_PICKING]";
	std::string HDR_ADV = "[ADVANCED]";
	std::string HDR_IDX = "[INDEXING]";
	std::string HDR_HLP = "[HELP]";
	std::string HDR_DEV = "[DEVELOPER]";

	// Type and Defaults are underlined
	// 'Required' is green
	// Options are in Bold
	// Name Type Required Descr Default
	//  20   12     10     47     20
	int name_w  = 18;
	int type_w  = 12;
	int req_w   = 9;
	int descr_w = 47; // description width
	int def_w = 20;
	ss << help_header << std::endl;
	for (auto opt : options)
	{
		std::string space_0 = "    "; // 4 chars
		std::string space_1 = "";
		std::string space_2 = "";

		req = std::get<3>(opt) ? "Required" : "Optional";
		pfx = std::get<0>(opt).size() > 1 ? "--" : "-";
		switch (std::get<2>(opt))
		{
		case OPT_CATEGORY::COMMON:
			if (std::get<3>(opt))
			{
				if (!set_req) {
					ss << space_0 << HDR_REQ << std::endl;
					set_req = true;
				}
			}
			else
			{
				if (!set_com)
				{
					ss << std::endl << space_0 << HDR_COM << std::endl;
					set_com = true;
				}
			}
			break;
		case OPT_CATEGORY::ADVANCED: 
			if (!set_adv) {
				ss << std::endl << space_0 << HDR_ADV << std::endl;
				set_adv = true;
			}
			break;
		case OPT_CATEGORY::OTU_PICKING:
			if (!set_otu) {
				ss << std::endl << space_0 << HDR_OTU << std::endl;
				set_otu = true;
			}
			break;
		case OPT_CATEGORY::HELP:
			if (!set_hlp) {
				ss << std::endl << space_0 << HDR_HLP << std::endl;
				set_hlp = true;
			}
			break;
		case OPT_CATEGORY::INDEXING:
			if (!set_idx) {
				ss << std::endl << space_0 << HDR_IDX << std::endl;
				set_idx = true;
			}
			break;
		case OPT_CATEGORY::DEVELOPER:
			if (!set_dev) {
				ss << std::endl << space_0 << HDR_DEV << std::endl;
				set_dev = true;
			}
			break;
		}

		int space_1_size = name_w - std::get<0>(opt).size() - pfx.size(); // name
		int space_2_size = type_w - std::get<1>(opt).size(); // type
		for (int cnt = 0; cnt < space_1_size; ++cnt)
			space_1 += " ";
		for (int cnt = 0; cnt < space_2_size; ++cnt)
			space_2 += " ";

		ss << space_0 << pfx << std::get<0>(opt) << space_1 << std::get<1>(opt) << space_2 << req << "  " << std::get<4>(opt) << std::endl;
	}
	std::cout << ss.str();
}

void Runopts::store_to_db(KeyValueDatabase &kvdb)
{}

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
		<< "  Copyright:    2016-2019 Clarity Genomics BVBA:" << std::endl
		<< "                Turnhoutseweg 30, 2340 Beerse, Belgium" << std::endl
		<< "                2014-2016 Knight Lab:" << std::endl
		<< "                Department of Pediatrics, UCSD, La Jolla" << std::endl
		<< "                2012-2014 Bonsai Bioinformatics Research Group:" << std::endl
		<< "                LIFL, University Lille 1, CNRS UMR 8022, INRIA Nord-Europe" << std::endl
		<< "  Disclaimer:   SortMeRNA comes with ABSOLUTELY NO WARRANTY; without even the" << std::endl
		<< "                implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE." << std::endl
		<< "                See the GNU Lesser General Public License for more details." << std::endl
		<< "  Contributors: Jenya Kopylova   jenya.kopylov@gmail.com" << std::endl
		<< "                Laurent Noé      laurent.noe@lifl.fr"      << std::endl
		<< "                Pierre Pericard  pierre.pericard@lifl.fr"  << std::endl
		<< "                Daniel McDonald  wasade@gmail.com"         << std::endl
		<< "                Mikaël Salson    mikael.salson@lifl.fr"    << std::endl
		<< "                Hélène Touzet    helene.touzet@lifl.fr"    << std::endl
		<< "                Rob Knight       robknight@ucsd.edu\n"     << std::endl;

	std::cout << ss.str();
}