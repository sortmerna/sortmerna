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
void help();
std::string get_user_home(); // util.cpp
unsigned int list_dir(std::string dpath);
bool dirExists(std::string dpath);
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
		ERR("option '--reads' requires a path to a reads FASTA/FASTQ file");
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
		std::cout << "File  [" << file << "] appears to be a relative path. Trying to locate in current directory and in work directory";
		fpath_a = std::filesystem::current_path() / file; // std::filesystem::path(file) makes no difference
		if (std::filesystem::exists(fpath_a))
		{
			std::cout << "File [" << fpath_a << "] exists" << std::endl;
		}
		else
		{
			fpath_a = std::filesystem::path(workdir) / file;
			if (!std::filesystem::exists(fpath_a))
			{
				std::cout << "File [" << fpath_a << "] does not exists" << std::endl;
			}
		}
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
		ERR("--ref must be followed by a file path (ex. --ref /path/to/file1.fasta)");
		exit(EXIT_FAILURE);
	}

	// verify the reference file exists and can be read
	// TODO:
	// if index already exists, no refs files are needed (only the names) =>
	// verify first the index exists based on the ref name. 
	// Verify physical presence and compare against metadata that can be stored in an index descriptor
	// or in the RocksDB DB (for this the options would need handle to the DB - OK?). 
	if (filesize(file) <= 0)
	{
		ss << STAMP << "File '" << file << "' either non-existent or empty or corrupt";
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}
	else
	{
		std::cout << STAMP << "'" << file << "'" << std::endl;
	}

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

	std::string basename = get_basename(file);
	// derive index file prefix from the reference file name
	// if we are here the Workdir is OK
	std::string idx_file_pfx = workdir + "/" + IDX_DIR + "/" + string_hash(basename);

	indexfiles.push_back(std::pair<std::string, std::string>(file, idx_file_pfx));
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
		std::cout << STAMP << "File name was not provided with option '--other [FILE]'. Using default name: [" << other_out_pfx << "]" << std::endl;
	}
	is_other = true;
} // ~Runopts::opt_other

void Runopts::opt_log(const std::string &val)
{
	if (write_log)
		WARN("'--" << OPT_LOG << "' is deprecated. True by default.");
} // ~Runopts::optLog

void Runopts::opt_de_novo_otu(const std::string &val)
{
	de_novo_otu = true;
} // ~Runopts::opt_de_novo_otu

void Runopts::opt_otu_map(const std::string &val)
{
	otumapout = true;
} // ~Runopts::opt_otu_map

void Runopts::opt_print_all_reads(const std::string &val)
{
	print_all_reads = true;
} // ~Runopts::optPrintAllReads

void Runopts::opt_pid(const std::string &val)
{
	pid = true;
} // ~Runopts::optPid

void Runopts::opt_paired_in(const std::string &val)
{
	if (pairedout)
	{
		ERR("'--paired_out' has been set, please choose one or the other, or use the default option");
		exit(EXIT_FAILURE);
	}

	pairedin = true;
} // ~Runopts::optPairedIn

void Runopts::opt_paired_out(const std::string &val)
{
	if (pairedin)
	{
		ERR("'--paired_in' has been set, please choose one or the other, or use the default option");
		exit(EXIT_FAILURE);
	}

	pairedout = true;
} // ~Runopts::optPairedOut

void Runopts::opt_match(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR("'--match [INT]' requires a positive integer as input (ex. --match 2)");
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
		ERR("--match [INT] has been set twice, please verify your choice");
		help();
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
		help();
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
		help();
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
		help();
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
		help();
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_num_seeds

/* --fastx */
void Runopts::opt_fastx(const std::string &val)
{
	if (is_fastxout)
	{
		ERR("--fastx has already been set once.");
		exit(EXIT_FAILURE);
	}
	else
	{
		is_fastxout = true;
	}
} // ~Runopts::opt_fastx

void Runopts::opt_sam(const std::string &val)
{
	if (samout)
	{
		ERR("--sam has already been set once.");
		exit(EXIT_FAILURE);
	}
	else
	{
		samout = true;
	}
} // ~Runopts::opt_sam

void Runopts::opt_blast(const std::string &val)
{
	std::stringstream ss;

	if (blastout)
	{
		ERR("--blast [STRING] has already been set once.");
		exit(EXIT_FAILURE);
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

	blastout = true;
} // ~Runopts::opt_blast

void Runopts::opt_min_lis(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR(": --min_lis [INT] requires an integer (>=0) as input (ex. --min_lis 2)."
			" Note: 0 signifies to search all high scoring reference sequences.");
		exit(EXIT_FAILURE);
	}

	// min_lis_gv has already been set
	if (min_lis_set)
	{
		ERR(": --min_lis [INT] has been set twice, please verify your choice.");
		help();
		exit(EXIT_FAILURE);
	}
	else
	{
		if ((sscanf(val.data(), "%d", &min_lis) != 1) || (min_lis < 0))
		{
			ERR(": --min_lis [INT] must be >= 0 (0 signifies to search all high scoring reference sequences).");
			exit(EXIT_FAILURE);
		}
		min_lis_set = true;
	}
} // ~Runopts::opt_min_lis

void Runopts::opt_best(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR(": --best [INT] requires an integer (> 0) as input (ex. --best 2).");
		exit(EXIT_FAILURE);
	}

	// best_set has already been set
	if (best_set)
	{
		ERR(" : --best [INT] has been set twice, please verify your choice.");
		help();
		exit(EXIT_FAILURE);
	}
	else
	{
		if ((sscanf(val.data(), "%d", &num_best_hits) != 1))
		{
			ERR(": could not read --best [INT] as integer");
			exit(EXIT_FAILURE);
		}
		best_set = true;
	}
} // ~Runopts::opt_best

void Runopts::opt_num_alignments(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR(": --num_alignments [INT] requires an integer (>=0) as input (ex. --num_alignments 2)"
			" Note: 0 signifies to output all alignments.");
		exit(EXIT_FAILURE);
	}

	if (num_alignments_set)
	{
		ERR(" :--num_alignments [INT] has been set twice, please verify your command parameters.");
		exit(EXIT_FAILURE);
	}

	// set number of alignments to output reaching the E-value
	num_alignments = atoi(val.data());
	if (num_alignments < 0)
	{
		ERR(": --num_alignments [INT] must be >= 0 (0 signifies to output all alignments).");
		exit(EXIT_FAILURE);
	}
	num_alignments_set = true;
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
		as_percent = true;

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
	full_search = true;
} // ~Runopts::opt_full_search

void Runopts::opt_SQ(const std::string &val)
{
	if (yes_SQ)
	{
		ERR(" : BOOL --SQ has been set twice, please verify your choice.");
		exit(EXIT_FAILURE);
	}

	yes_SQ = true;
} // ~Runopts::optSQ

void Runopts::opt_passes(const std::string &val)
{
	if (passes_set)
	{
		ERR(" : --passes [INT,INT,INT] has been set twice, please verify your choice.");
		exit(EXIT_FAILURE);
	}

	// set passes
	for (auto pos = val.find(",", 0), pos_from = pos-pos, count = pos-pos; pos != std::string::npos; pos = val.find(",", pos_from))
	{
		auto tok = val.substr(pos_from, pos);
		pos_from += pos +1;
		if (++count > 3) 
		{
			ERR(" : exactly 3 integers has to be provided with '--passes [INT,INT,INT]'");
			exit(EXIT_FAILURE);
		}
		auto skiplen = std::stoi(tok);
		if (skiplen > 0)
			skiplengths.emplace_back(skiplen);
		else
		{
			ERR(" : all three integers in --passes [INT,INT,INT] "
				"must contain positive integers where 0<INT<(shortest read length).");
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
	help();
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
	if (!forward)
	{
		forward = true;
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
	if (!reverse)
	{
		reverse = true;
	}
	else
	{
		ERR(" : BOOL '-R' has been set more than once, please check your command parameters.");
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_R

/* Help */
void Runopts::opt_h(const std::string &val)
{
	about();
	help();
	exit(0);
} // ~Runopts::opt_h

void Runopts::opt_v(const std::string &val)
{
	verbose = true;
} // ~Runopts::opt_v

void Runopts::opt_N(const std::string &val)
{
	// match ambiguous N's
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
	if (val.size() == 0)
	{
		ERR(": -a [INT] requires an integer for number of Processor threads (ex. -a 8)");
		exit(EXIT_FAILURE);
	}

	num_proc_thread = std::stoi(val);
} // ~Runopts::opt_a_numProcThreads

  /* Number of threads to use */
void Runopts::opt_threads(const std::string &val)
{
	if (val.size() == 0)
	{
		ERR(": --threads [INT:INT:INT] requires 3 integers for number of "
			<< "Read:Write:Processor threads (ex. --threads 1:1:8)");
		exit(EXIT_FAILURE);
	}

	std::istringstream strm(val);
	std::string tok;
	for (int i = 0; std::getline(strm, tok, ':'); ++i)
	{
		switch (i)
		{
		case 0: num_read_thread = std::stoi(tok); break;
		case 1: num_write_thread = std::stoi(tok); break;
		case 2: num_proc_thread = std::stoi(tok); break;
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
	dbg_put_kvdb = true;
}

void Runopts::opt_default(const std::string &opt)
{
	std::stringstream ss;
	ss << STAMP << "Option: '" << opt << "' is not recognized";
	ERR(ss.str());
	help();
	exit(EXIT_FAILURE);
} // ~Runopts::opt_default

  /* Processing task */
void Runopts::opt_task(const std::string &val)
{
	int taskOpt = std::stoi(val);

	if (taskOpt > 4) 
	{
		std::stringstream ss;
		ss << "Option '−−task' can only take values in range [0..4] Provided value is " << taskOpt;
		ERR(ss.str());
		exit(EXIT_FAILURE);
	}

	switch (taskOpt)
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
	interactive = true;
} // ~Runopts::optInteractive

/* Work directory setup */
void Runopts::opt_workdir(const std::string &path)
{
	if (path.size() == 0)
	{
		workdir = get_user_home() + "/sortmerna";
		std::cout << "'workdir' option not provided. Using USERDIR to set the working directory: [" << workdir << "]" << std::endl;
	}
	else
		workdir = path;
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
		ss << "Option 'L' given without value. Using default: " << lnwin_gv;
		WARN(ss.str());
	}
	else
	{
		int lnwin_t = std::stoi(val);

		if (lnwin_t <= 0 || lnwin_t % 2 == 1 || lnwin_t < 8 || lnwin_t > 26)
		{
			ss.str("");
			ss << STAMP 
				<< "Option L takes a Positive Even integer between 8 and 26 inclusive e.g. 10, 12, 14, .. , 20. Provided value: " 
				<< lnwin_t << " Default will be used: " << lnwin_gv;
			WARN(ss.str());
		}
		else
		{
			lnwin_gv = lnwin_t;
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

	std::cout << STAMP << "Key-value DB location (" << kvdbPath << ")" << std::endl;

	if (dirExists(kvdbPath))
	{
		// dir exists and not empty
		auto count = list_dir(kvdbPath);
		if (count > 0) 
		{
			// TODO: Store some metadata in DB to verify the alignment.
			// kvdb.verify()
			if (ALIGN_REPORT::align == alirep || ALIGN_REPORT::all == alirep || ALIGN_REPORT::alipost == alirep)
			{
				// if (kvdb.verify()) // TODO
				std::cout << STAMP << "Database (" << kvdbPath << ") exists. Alignment information OK" << std::endl;
			}
		}
		else
		{
			if (ALIGN_REPORT::postproc == alirep || ALIGN_REPORT::report == alirep)
			{
				std::cout << STAMP << "KVDB (" << kvdbPath << " is empty. Will run alignment" << std::endl;
			}
			// dir exists and empty -> use
		}
	}
	else
	{
		// dir does not exist -> try creating
		std::cout << STAMP << "Database (" << kvdbPath << ") will be created" << std::endl;
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
		verbose = true;
		about();
		ERR("Missing required command options");
		exit(EXIT_FAILURE);
	}

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

	// check required options were provided
	for (auto opt : options)
	{
		if (std::get<0>(opt.second))
		{
			if (mopt.count(opt.first) == 0)
			{
				std::cout << "Missing required flag: " << opt.first << std::endl;
			}
		}
	}

	// Process options
	// process WORKDIR first as other options depend on it
	auto wd_it = mopt.find("workdir");
	std::string wdir = "";
	if (wd_it != mopt.end())
	{
		wdir = wd_it->second;
		mopt.erase(wd_it); // remove to prevent repeated processing below
	}
	opt_workdir(wdir);

	// loop through the rest of options
	for (auto opt : mopt)
	{
		std::cout << STAMP << "Processing option: " << opt.first << " with value: " << opt.second << std::endl;

		int count = options.count(opt.first);
		if (count > 0) 
		{
			//std::string descr = std::get<1>(options.at(opt.first));
			//std::get<2>(options.at(opt.first))(opt.second); // call processing function for the given option
			std::invoke(std::get<2>(options.at(opt.first)), this, opt.second);
		}
		// passed option is not recognized
		else
		{
			opt_default(opt.first);
		}
	}

	// validate the options. TODO: should be part of options validation
	test_kvdb_path();
	validate();
} // ~Runopts::process

/**
 * Validate the options and setup defaults
 */
void Runopts::validate()
{
	// No output format has been chosen
	if (!(is_fastxout || blastout || samout || otumapout || de_novo_otu))
	{
		blastout = true;
		std::cout << STAMP << "No output format has been chosen (fastx/sam/blast/otu_map/log), Using default blast" << std::endl;
	}

	// Options --paired_in and --paired_out can only be used with FASTA/Q output
	if (!is_fastxout && (pairedin || pairedout))
	{
		ERR(": options '--paired_in' and '--paired_out' must be accompanied by option '--fastx'.\n"
			"  These BOOLs are for FASTA and FASTQ output files to maintain paired reads together.");
		exit(EXIT_FAILURE);
	}

	// An OTU map can only be constructed with the single best alignment per read
	if (otumapout && num_alignments_set)
	{
		ERR("'--otu_map' cannot be set together with --num_alignments [INT].\n"
			"   The option --num_alignments [INT] doesn't keep track of"
			" the best alignment which is required for constructing an OTU map.\n"
			"   Use --otu_map with --best [INT] instead.");
		exit(EXIT_FAILURE);
	}

	// If --num_alignments output was chosen, check an alignment format has also been chosen
	if (num_alignments_set && !(blastout || samout || is_fastxout))
	{
		ERR(" : --num_alignments [INT] has been set but no output "
			"format has been chosen (--blast, --sam or --fastx).");
		exit(EXIT_FAILURE);
	}

	// If --best output was chosen, check an alignment format has also been chosen
	if (best_set && !(blastout || samout || otumapout))
	{
		ERR(" : --best [INT] has been set but no output "
			"format has been chosen (--blast or --sam or --otu_map).");
		exit(EXIT_FAILURE);
	}

	// Check gap extend score < gap open score
	if (gap_extension > gap_open)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --gap_ext [INT] must be less than --gap_open [INT].\n\n",
			RED, COLOFF, __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	// Option --print_all_reads can only be used with Blast-like tabular
	// and SAM formats (not pairwise)
	if (print_all_reads && blastout && blastFormat != BlastFormat::TABULAR)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --print_all_reads [BOOL] can only be used with BLAST-like "
			"tabular format.\n\n", RED, COLOFF, __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	// Only one of these options is allowed (--best outputs one alignment,
	// --num_alignments outputs > 1 alignments)
	if (best_set && num_alignments_set)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --best [INT] and --num_alignments [INT] cannot "
			"be set together. \n", RED, COLOFF, __LINE__, __FILE__);
		fprintf(stderr, "  (--best [INT] will search INT highest scoring reference sequences "
			"and output a single best alignment, whereas --num_alignments [INT] will "
			"output the first INT alignments).\n\n");
	}

	// Option --min_lis [INT] can only accompany --best [INT]
	if (min_lis_set && num_alignments_set)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --min_lis [INT] and --num_alignments [INT] cannot "
			"be set together. \n", RED, COLOFF, __LINE__, __FILE__);
		fprintf(stderr, "  --min_lis [INT] can only be used with --best [INT] (refer to "
			"the User manual on this option).\n\n");
		exit(EXIT_FAILURE);
	}

	// Option --mis_lis INT accompanies --best INT, cannot be set alone
	if (min_lis_set && !best_set)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --min_lis [INT] must be set together with --best "
			"[INT].\n\n", RED, COLOFF, __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	// %id and %coverage can only be used with --otu_map
	if (((align_id > 0) || (align_cov > 0)) && !otumapout)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --id [INT] and --coverage [INT] can only be used "
			"together with --otu_map.\n", RED, COLOFF, __LINE__, __FILE__);
		fprintf(stderr, "  These two options are used for constructing the OTU map by "
			"filtering alignments passing the E-value threshold.\n\n");
		exit(EXIT_FAILURE);
	}

	// the list of arguments is correct, welcome the user!
	if (verbose) about();
	// if neither strand was selected for search, search both
	if (!forward && !reverse)
	{
		forward = true;
		reverse = true;
	}
	// default number of threads is 1
	//if (numcpu_gv < 0) numcpu_gv = 1;
	// default E-value
	if (evalue < 0.0) evalue = 1;
	// SW alignment parameters
	if (!match_set) match = 2;
	if (!mismatch_set) mismatch = -3;
	if (!gap_open_set) gap_open = 5;
	if (!gap_ext_set) gap_extension = 2;
	if (!match_ambiguous_N) score_N = mismatch;

	// default method for searching alignments
	if (!best_set && !num_alignments_set)
	{
		// FASTA/FASTQ output, stop searching for alignments after the first match
		if (is_fastxout && !(blastout || samout || otumapout || write_log || de_novo_otu))
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
	if (best_set && !min_lis_set) min_lis = 2;
	// default number of seed hits before searching for candidate LIS
	if (seed_hits < 0) seed_hits = 2;
	// default number of nucleotides to add to each edge of an alignment
	// region before extension
	if (edges < 0) edges = 4;
	// activate heuristic for stopping search (of 1-error matches) after
	// finding 0-error match
	if (!full_search_set) full_search = false;
	// default %id to keep alignment
	if (align_id < 0)
	{
		// if OTU-map is chosen, set default similarity to 0.97
		if (otumapout) align_id = 0.97;
		else align_id = 0;
	}
	// default %query coverage to keep alignment
	if (align_cov < 0)
	{
		// if OTU-map is chosen, set default coverage to 0.97
		if (otumapout) align_cov = 0.97;
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
		<< "  Contributors: Jenya Kopylova   jenya.kopylov@gmail.com " << std::endl
		<< "                Laurent Noé      laurent.noe@lifl.fr"      << std::endl
		<< "                Pierre Pericard  pierre.pericard@lifl.fr"  << std::endl
		<< "                Daniel McDonald  wasade@gmail.com"         << std::endl
		<< "                Mikaël Salson    mikael.salson@lifl.fr"    << std::endl
		<< "                Hélène Touzet    helene.touzet@lifl.fr"    << std::endl
		<< "                Rob Knight       robknight@ucsd.edu\n"     << std::endl;

	std::cout << ss.str();
}



/*! @fn printlist()
*  @brief outputs options menu
*  @param none
*  @return none
*/
void help()
{
	std::stringstream ss;

	ss << std::endl
		<< "  usage:   sortmerna --ref REFERENCE_FILE --reads READS_FILE [OPTIONS]:" << std::endl
		<< std::endl
		<< "  -------------------------------------------------------------------------------------------------------------"  << std::endl
		<< "  | option              type-format       description                                              default    |"  << std::endl
		<< "  -------------------------------------------------------------------------------------------------------------"  << std::endl
		<< "  [REQUIRED OPTIONS]: "                                                                                           << std::endl << BOLD
		<< "    --ref            "                                                                                            << COLOFF << UNDL
		<<                       "  STRING       "                                                                            << COLOFF
		<<                                       "   FASTA reference file                                      "              << GREEN 
		<<                                                                                                     "mandatory"    << COLOFF << std::endl
		<< "                                         Use multiple 'ref' options to specify multiple            "              << std::endl
		<< "                                         reference files e.g. --ref FILE_1 --ref FILE_2 ...        "              << std::endl << BOLD
		<< "    --reads          "                                                                                            << COLOFF << UNDL
		<<                       "  STRING       "                                                                            << COLOFF
		<<                                       "   FASTA/FASTQ raw reads file                                "              << GREEN 
		<<                                                                                                     "mandatory"    << COLOFF << std::endl
		<< "    --aligned        "                                                                                            << COLOFF << UNDL
		<<                       "  STRING       "                                                                            << COLOFF
		<<                                       "   aligned reads filepath + base file name                   "              << GREEN 
		<<                                                                                                     "mandatory"    << COLOFF << std::endl
		<< "                                         (appropriate extension will be added)"                                   << std::endl << std::endl
		<< "  [COMMON OPTIONS]: "                                                                                             << std::endl << BOLD
		<< "    --other         "                                                                                             << COLOFF << UNDL
		<<                      "  STRING        "                                                                            << COLOFF
		<<                                       "   rejected reads filepath + base file name                  "              << std::endl
		<< "                                         (appropriate extension will be added)                     "              << std::endl << BOLD
		<< "    --fastx         "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   output FASTA/FASTQ file                                   "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
		<< "                                         (for aligned and/or rejected reads)"                                     << std::endl << BOLD
		<< "    --sam           "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   output SAM alignment                                      "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
		<< "                                         (for aligned reads only)"                                                << std::endl << BOLD
		<< "    --SQ            "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   add SQ tags to the SAM file                               "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl << BOLD
		<< "    --blast         "                                                                                             << COLOFF << UNDL
		<<                      "  STRING        "                                                                            << COLOFF
		<<                                       "   output alignments in various Blast-like formats"                         << std::endl
		<< "                                           '0' - pairwise"                                                        << std::endl
		<< "                                           '1' - tabular (Blast -m 8 format)"                                     << std::endl
		<< "                                           '1 cigar' - tabular + column for CIGAR "                               << std::endl
		<< "                                           '1 cigar qcov' - tabular + columns for CIGAR"                          << std::endl
		<< "                                                      and query coverage"                                         << std::endl
		<< "                                           '1 cigar qcov qstrand' - tabular + columns for CIGAR,"                 << std::endl
		<< "                                                                query coverage and strand"                        << std::endl << BOLD
		<< "    --log           "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   output overall statistics                                 "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
#ifdef NOMASK_option
		<< "     " << BOLD << "--no-mask" << COLOFF << "         " << UNDL << "BOOL" << COLOFF
		<< "            do not mask low occurrence (L/2)-mers when searching           " << UNDL << "on" << COLOFF << std::endl
		<< "                                       for seeds of length L\n" << std::endl
#endif
		                                                                                                                      << BOLD
		<< "    --num_alignments"                                                                                             << COLOFF << UNDL
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   report first INT alignments per read reaching E-value     "              << UNDL 
		<<                                                                                                     "-1"           << COLOFF << std::endl
		<< "                                        (--num_alignments 0 signifies all alignments will be output)"             << std::endl << COLOFF
		<< "       OR"                                                                                                        << GREEN
		<<           " (default)"                                                                                             << std::endl << COLOFF << BOLD
		<< "    --best          "                                                                                             << COLOFF << UNDL
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   report INT best alignments per read reaching E-value      "              << UNDL 
		<<                                                                                                     "1"            << COLOFF << std::endl
		<< "                                         by searching --min_lis INT candidate alignments"                         << std::endl
		<< "                                        (--best 0 signifies all candidate alignments will be searched)"           << std::endl << BOLD
		<< "    --min_lis       "                                                                                             << COLOFF << UNDL
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   search all alignments having the first INT longest LIS    "              << UNDL 
		<<                                                                                                     "2"            << COLOFF << std::endl
		<< "                                         LIS stands for Longest Increasing Subsequence, it is "                   << std::endl
		<< "                                         computed using seeds' positions to expand hits into"                     << std::endl
		<< "                                         longer matches prior to Smith-Waterman alignment. "                      << std::endl << BOLD
		<< "   --print_all_reads"                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   output null alignment strings for non-aligned reads       "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
		<< "                                         to SAM and/or BLAST tabular files"                                       << std::endl << BOLD
		<< "    --paired_in     "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   both paired-end reads go in --aligned fasta/q file        "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
		<< "                                        (interleaved reads only, see Section 4.2.4 of User Manual)"               << std::endl << BOLD
		<< "    --paired_out    "                                                                                             << COLOFF << UNDL
		<<                      "   BOOL         "                                                                            << COLOFF
		<<                                       "   both paired-end reads go in --other fasta/q file          "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
		<< "                                       (interleaved reads only, see Section 4.2.4 of User Manual)"                << std::endl << BOLD
		<< "    --match         "                                                                                             << COLOFF << UNDL
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   SW score (positive integer) for a match                   "              << UNDL 
		<<                                                                                                     "2"            << COLOFF << std::endl << BOLD
		<< "    --mismatch      "                                                                                             << COLOFF << UNDL 
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   SW penalty (negative integer) for a mismatch              "              << UNDL 
		<<                                                                                                     "-3"           << COLOFF << std::endl << BOLD
		<< "    --gap_open      "                                                                                             << COLOFF << UNDL 
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   SW penalty (positive integer) for introducing a gap       "              << UNDL 
		<<                                                                                                     "5"            << COLOFF << std::endl << BOLD
		<< "    --gap_ext       "                                                                                             << COLOFF << UNDL
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   SW penalty (positive integer) for extending a gap         "              << UNDL 
		<<                                                                                                     "2"            << COLOFF << std::endl << BOLD
		<< "    -N              "                                                                                             << COLOFF << UNDL 
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   SW penalty for ambiguous letters (N's)                    "              << std::endl
		<< "                                         scored as --mismatch"                                                    << std::endl << BOLD
		<< "    -F              "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   search only the forward strand                            "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl << BOLD
		<< "    -R              "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   search only the reverse-complementary strand              "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl << BOLD
		<< "    -e              "                                                                                             << COLOFF << UNDL
		<<                      "  DOUBLE        "                                                                            << COLOFF
		<<                                       "   E-value threshold                                         "              << UNDL 
		<<                                                                                                     "1"            << COLOFF << std::endl << BOLD
		<< "    -v              "                                                                                             << COLOFF << UNDL 
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   verbose                                                   "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl << std::endl
		<< "  [OTU PICKING OPTIONS]: "                                                                                        << std::endl << BOLD
		<< "    --id            "                                                                                             << COLOFF << UNDL 
		<<                      "  DOUBLE        "                                                                            << COLOFF
		<<                                       "   %%id similarity threshold (the alignment must             "              << UNDL 
		<<                                                                                                     "0.97"         << COLOFF << std::endl
		<< "                                         still pass the E-value threshold)\n"                                     << std::endl << BOLD
		<< "    --coverage      "                                                                                             << COLOFF << UNDL
		<<                      "  DOUBLE        "                                                                            << COLOFF
		<<                                       "   %%query coverage threshold (the alignment must            "              << UNDL 
		<<                                                                                                     "0.97"         << COLOFF << std::endl
		<< "                                         still pass the E-value threshold)\n"                                     << std::endl << BOLD
		<< "    --de_novo_otu   "                                                                                             << COLOFF << UNDL 
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   FASTA/FASTQ file for reads matching database < %%id       "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
		<< "                                         (set using --id) and < %%cov (set using --coverage) "                    << std::endl
		<< "                                         (alignment must still pass the E-value threshold)"                       << std::endl << BOLD
		<< "    --otu_map       "                                                                                             << COLOFF << UNDL 
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   output OTU map (input to QIIME's make_otu_table.py)       "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl << std::endl
		<< "  [ADVANCED OPTIONS] (see SortMeRNA user manual for more details): "                                              << std::endl << BOLD
		<< "    --passes        "                                                                                             << COLOFF << UNDL 
		<<                      "  INT,INT,INT   "                                                                            << COLOFF
		<<                                       "   three intervals at which to place the seed on the read    "              << UNDL 
		<<                                                                                                     "L,L/2,3"      << COLOFF << std::endl
		<< "                                         (L is the seed length set in ./indexdb_rna)"                             << std::endl << BOLD
		<< "    --edges         "                                                                                             << COLOFF << UNDL 
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   number (or percent if INT followed by %% sign) of         "              << UNDL 
		<<                                                                                                     "4"            << COLOFF << std::endl
		<< "                                         nucleotides to add to each edge of the read"                             << std::endl
		<< "                                         prior to SW local alignment "                                            << std::endl << BOLD
		<< "    --num_seeds     "                                                                                             << COLOFF << UNDL 
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   number of seeds matched before searching                  "              << UNDL 
		<<                                                                                                     "2"            << COLOFF << std::endl
		<< "                                         for candidate LIS "                                                      << std::endl << BOLD
		<< "    --full_search   "                                                                                             << COLOFF << UNDL 
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   search for all 0-error and 1-error seed                   "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl
		<< "                                         matches in the index rather than stopping"                               << std::endl
		<< "                                         after finding a 0-error match (<1%% gain in"                             << std::endl
		<< "                                         sensitivity with up four-fold decrease in speed)"                        << std::endl << BOLD
		<< "    --pid           "                                                                                             << COLOFF << UNDL 
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   add pid to output file names                              "              << UNDL 
		<<                                                                                                     "off"          << COLOFF << std::endl << std::endl
		<< "   [HELP]: "                                                                                                      << std::endl << BOLD
		<< "    -h              "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF 
		<<                                       "   help"                                                                    << std::endl << BOLD
		<< "    --version       "                                                                                             << COLOFF << UNDL 
		<<                      "  BOOL          "                                                                            << COLOFF 
		<<                                       "   SortMeRNA version number"                                                << std::endl << std::endl
		<< "  [DEVELOPER OPTIONS]: "                                                                                          << std::endl << BOLD
		<< "    --cmd           "                                                                                             << COLOFF << UNDL
		<<                      "  BOOL          "                                                                            << COLOFF
		<<                                       "   launch an interactive session (command prompt)"                          << std::endl << BOLD
		<< "    --task          "                                                                                             << COLOFF << UNDL
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   Processing Task"                                                         << std::endl
		<< "                                           0 - align Only perform alignment"                                      << std::endl
		<< "                                           1 - post-processing (log writing)"                                     << std::endl
		<< "                                           2 - generate reports"                                                  << std::endl
		<< "                                           3 - align and post−process"                                            << std::endl
		<< "                                           4 - all (default)"                                                     << std::endl << std::endl << BOLD
		<< "    -d              "                                                                                             << COLOFF << UNDL
		<<                      "  STRING        "                                                                            << COLOFF
		<<                                       "   key-value datastore FULL folder path                      "              << UNDL
		<<                                                                                                     "USERDIR/kvdb" << COLOFF << std::endl << BOLD
		<< "    -a              "                                                                                             << COLOFF << UNDL
		<<                      "  INT           "                                                                            << COLOFF
		<<                                       "   number of threads to use                                  "              << UNDL
		<<                                                                                                     "numCores"     << COLOFF << std::endl << BOLD
		<< "    --threads       "                                                                                             << COLOFF << UNDL
		<<                      "  INT:INT:INT   "                                                                            << COLOFF
		<<                                       "   number of Read:Write:Process threads to use               "              << UNDL
		<<                                                                                                     "1:1:numCores" << COLOFF << std::endl << BOLD
		<< "    --thpp          "                                                                                             << COLOFF << UNDL
		<<                      "  INT:INT:INT   "                                                                            << COLOFF
		<<                                       "   number of Post-Processing Read:Process threads to use     "              << UNDL
		<<                                                                                                     "1:1"          << COLOFF << std::endl << BOLD
		<< "    --threp         "                                                                                             << COLOFF << UNDL
		<<                      "  INT:INT:INT   "                                                                            << COLOFF
		<<                                       "   number of Report Read:Process threads to use              "              << UNDL
		<<                                                                                                     "1:1"          << COLOFF << std::endl << std::endl;
		
	std::cout << ss.str();
}//~help()