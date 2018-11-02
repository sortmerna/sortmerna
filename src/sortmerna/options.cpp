/* 
 * FILE: options.cpp
 * Created: Jun 07, 2018 Thu
 */

 // TODO: BUG: if SMR headers moved down after 3rd party, 'timeval' struct gets 'redefined' - compiler error. That's a header mess bug.

 // 3rd party
#include "zlib.h"

#include "version.h"
#include "build_version.h"
#include "options.hpp"
#include "common.hpp"

 // standard
#include <limits>
#include <dirent.h>
#include <unistd.h>
#include <sstream>
#include <iostream>
#include <cstring> // strerror, strrchr, memcpy, strcpy, strpbrk
#include <fcntl.h>


#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

/*! @brief Measure time using this variable. */
//timeval t;

// forward
void welcome();
void printlist();
std::string get_user_home();
unsigned int list_dir(std::string dpath);
bool dirExists(std::string dpath);

void Runopts::optReads(char **argv, int &narg)
{
	if (have_reads_gz)
	{
		fprintf(stderr, "\n %sERROR%s: option --reads-gz has also been set, only one of "
			"--reads-gz or --reads is permitted\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: a path to a reads FASTA/FASTQ file "
			"must be given after the option --reads\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		// check the file exists
		if (FILE *file = fopen(argv[narg + 1], "rb"))
		{
			// get size of file
			fseek(file, 0, SEEK_END);
			size_t filesize = ftell(file);

			// set exit BOOL to exit program after outputting
			// empty files, sortmerna will not execute after
			// that call (in paralleltraversal.cpp)
			if (!filesize) exit_early = true;
			// reset file pointer to start of file
			fseek(file, 0, SEEK_SET);

			readsfile = argv[narg + 1];
			narg += 2;
			fclose(file);

			have_reads = true;
		}
		else
		{
			fprintf(stderr, "\n  %sERROR%s: the file %s could not be opened: "
				"%s.\n\n", RED, COLOFF, argv[narg + 1], strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
} // ~Runopts::optReads

void Runopts::optReadsGz(char **argv, int &narg)
{
	if (have_reads)
	{
		fprintf(stderr, "\n %sERROR%s: option --reads has also been set, only one of "
			"--reads or --reads-gz is permitted\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: a path to a reads FASTA/FASTQ compressed (.zip, .gz) file "
			"must be given after the option --reads-gz\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		// check the file exists
		if (gzFile file = gzopen(argv[narg + 1], "r"))
		{
			// get size of file
			gzseek(file, 0, SEEK_END);
			size_t filesize = gztell(file);

			// set exit BOOL to exit program after outputting
			// empty files, sortmerna will not execute after
			// that call (in paralleltraversal.cpp)
			if (!filesize) exit_early = true;
			// reset file pointer to start of file
			gzseek(file, 0, SEEK_SET);

			readsfile = argv[narg + 1];
			narg += 2;
			gzclose(file);

			have_reads_gz = true;
		}
		else
		{
			fprintf(stderr, "\n  %sERROR%s: the file %s could not be opened: "
				"%s.\n\n", RED, COLOFF, argv[narg + 1], strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
} // ~Runopts::optReadsGz

void Runopts::optRef(char **argv, int &narg)
{
	std::stringstream ss;

	if (argv[narg + 1] == NULL)
	{
		ss << "\n  " << RED << "ERROR" << COLOFF
			<< ": --ref must be followed by at least one entry (ex. --ref /path/to/file1.fasta,/path/to/index1)"
			<< std::endl << std::endl;
		std::cerr << ss.str(); ss.str("");
		exit(EXIT_FAILURE);
	}

	// check path
	char *ptr = argv[narg + 1];
	while (*ptr != '\0')
	{
		// get the FASTA file path + name
		char fastafile[2000];
		char *ptr_fastafile = fastafile;

		// the reference database FASTA file
		while (*ptr != ',' && *ptr != '\0')
		{
			*ptr_fastafile++ = *ptr++;
		}
		*ptr_fastafile = '\0';
		if (*ptr == '\0')
		{
			fprintf(stderr, "   %sERROR%s: the FASTA reference file name %s must be followed "
				" by an index name.\n\n", RED, COLOFF, fastafile);
			exit(EXIT_FAILURE);
		}
		ptr++; //skip the ',' delimiter

			   // check reference FASTA file exists & is not empty
		if (FILE *file = fopen(fastafile, "rb"))
		{
			// get file size
			fseek(file, 0, SEEK_END);
			size_t filesize = ftell(file);
			if (!filesize) exit_early = true;
			// reset file pointer to start of file
			fseek(file, 0, SEEK_SET);
			fclose(file);
		}
		else
		{
			fprintf(stderr, "\n  %sERROR%s: the file %s could not be opened: "
				" %s.\n\n", RED, COLOFF, fastafile, strerror(errno));
			exit(EXIT_FAILURE);
		}

		// get the index path + name
		char indexfile[2000];
		char *ptr_indexfile = indexfile;
		// the reference database index name
		while (*ptr != DELIM && *ptr != '\0') *ptr_indexfile++ = *ptr++;
		*ptr_indexfile = '\0';
		if (*ptr != '\0') ptr++; //skip the ':' delimiter

								 // check the directory where to write the index exists
		char dir[500];
		char *ptr_end = strrchr(indexfile, '/');
		if (ptr_end != NULL)
		{
			memcpy(dir, indexfile, (ptr_end - indexfile));
			dir[(int)(ptr_end - indexfile)] = '\0';
		}
		else
		{
			strcpy(dir, "./");
		}

		if (DIR *dir_p = opendir(dir)) closedir(dir_p);
		else
		{
			if (ptr_end != NULL)
				fprintf(stderr, "\n  %sERROR%s: the directory %s for writing index "
					"'%s' could not be opened. The full directory path must be "
					"provided (ex. no '~'). \n\n", RED, COLOFF,
					dir, ptr_end + 1);
			else
				fprintf(stderr, "\n  %sERROR%s: the directory %s for writing index "
					"'%s' could not be opened. The full directory path must be "
					"provided (ex. no '~'). \n\n", RED, COLOFF,
					dir, indexfile);

			exit(EXIT_FAILURE);
		}

		// check index file names are distinct
		for (int i = 0; i < (int)indexfiles.size(); i++)
		{
			if ((indexfiles[i].first).compare(fastafile) == 0)
			{
				fprintf(stderr, "\n  %sWARNING%s: the FASTA file %s has been entered "
					"twice in the list. It will be searched twice. "
					"\n\n", "\033[0;33m", COLOFF, fastafile);
			}
			else if ((indexfiles[i].second).compare(indexfile) == 0)
			{
				fprintf(stderr, "\n  %sWARNING%s: the index name %s has been entered "
					"twice in the list. It will be searched twice.\n\n", "\033[0;33m",
					COLOFF, indexfile);
			}
		}

		indexfiles.push_back(std::pair<std::string, std::string>(fastafile, indexfile));

	}//~while (*ptr != '\0')

	narg += 2;
} // ~Runopts::optRef

void Runopts::optAligned(char **argv, int &narg)
{
	if ((argv[narg + 1] == NULL) || (argv[narg + 1][0] == '-'))
	{
		fprintf(stderr, "\n  %sERROR%s: a filename must follow the option --aligned [STRING]\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		// check if the directory where to write exists
		char dir[500];
		char *ptr = strrchr(argv[narg + 1], '/');
		if (ptr != NULL)
		{
			memcpy(dir, argv[narg + 1], (ptr - argv[narg + 1]));
			dir[(int)(ptr - argv[narg + 1])] = '\0';
		}
		else
		{
			strcpy(dir, "./");
		}

		if (DIR *dir_p = opendir(dir))
		{
			filetype_ar.assign(argv[narg + 1]);
			narg += 2;
			closedir(dir_p);
		}
		else
		{
			fprintf(stderr, "\n  %sERROR%s: the --aligned <STRING> directory "
				"%s could not be opened: %s.\n\n", RED, COLOFF,
				dir, strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
} // ~Runopts::optAligned

void Runopts::optOther(char **argv, int &narg)
{
	if ((argv[narg + 1] == NULL) || (argv[narg + 1][0] == '-'))
	{
		fprintf(stderr, "\n  %sERROR%s: a filename must follow the option "
			"--other [STRING]\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		// check if the directory where to write exists
		char dir[500];
		char *ptr = strrchr(argv[narg + 1], '/');
		if (ptr != NULL)
		{
			memcpy(dir, argv[narg + 1], (ptr - argv[narg + 1]));
			dir[(int)(ptr - argv[narg + 1])] = '\0';
		}
		else
		{
			strcpy(dir, "./");
		}

		if (DIR *dir_p = opendir(dir))
		{
			filetype_or.assign(argv[narg + 1]);
			narg += 2;
			closedir(dir_p);
		}
		else
		{
			fprintf(stderr, "\n  %sERROR%s: the --other %s directory could not be opened, please check it exists.\n\n", RED, COLOFF, dir);
			exit(EXIT_FAILURE);
		}
	}
} // ~Runopts::optOther

void Runopts::optLog(char **argv, int &narg)
{
	if (doLog)
	{
		fprintf(stderr, "\n  %sERROR%s: --log has already been set once.\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		doLog = true;
		narg++;
	}
} // ~Runopts::optLog

void Runopts::optDeNovoOtu(char **argv, int &narg)
{
	if (de_novo_otu)
	{
		fprintf(stderr, "\n  %sERROR%s: --de_novo_otu has already been set once.\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		de_novo_otu = true;
		narg++;
	}
} // ~Runopts::optDeNovoOtu

void Runopts::optOtuMap(char **argv, int &narg)
{
	if (otumapout)
	{
		fprintf(stderr, "\n  %sERROR%s: --otu_map has already been set once.\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		otumapout = true;
		narg++;
	}
} // ~Runopts::optOtuMap

void Runopts::optPrintAllReads(char **argv, int &narg)
{
	if (print_all_reads)
	{
		fprintf(stderr, "\n  %sERROR%s: --print_all_reads has already been set once.\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		print_all_reads = true;
		narg++;
	}
} // ~Runopts::optPrintAllReads

void Runopts::optPid(char **argv, int &narg)
{
	if (pid)
	{
		fprintf(stderr, "\n  %sERROR%s: --pid has already been set once.\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		pid = true;
		narg++;
	}
} // ~Runopts::optPid

void Runopts::optPairedIn(char **argv, int &narg)
{
	if (pairedin)
	{
		fprintf(stderr, "\n  %sERROR%s: --paired_in has already been set once.\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else if (pairedout)
	{
		fprintf(stderr, "\n  %sERROR%s: --paired_out has been set, please choose "
			"one or the other, or use the default option.\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		pairedin = true;
		narg++;
	}
} // ~Runopts::optPairedIn

void Runopts::optPairedOut(char **argv, int &narg)
{
	if (pairedout)
	{
		fprintf(stderr, "\n  %sERROR%s: --paired_out has already been set once.\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else if (pairedin)
	{
		fprintf(stderr, "\n %sERROR%s: --paired_in has been set, please choose one "
			"or the other, or use the default option.\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		pairedout = true;
		narg++;
	}
} // ~Runopts::optPairedOut

void Runopts::optMatch(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --match [INT] requires a positive integer as "
			"input (ex. --match 2).\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// set match
	if (!match_set)
	{
		match = atoi(argv[narg + 1]);
		narg += 2;
		match_set = true;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: --match [INT] has been set twice, please "
			"verify your choice\n\n", RED, COLOFF);
		printlist();
	}
} // ~Runopts::optMatch

void Runopts::optMismatch(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --mismatch [INT] requires a negative integer "
			"input (ex. --mismatch -2)\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// set mismatch
	if (!mismatch_set)
	{
		mismatch = atoi(argv[narg + 1]);
		if (mismatch > 0)
		{
			fprintf(stderr, "\n  %sERROR%s: --mismatch [INT] requires a negative "
				"integer input (ex. --mismatch -2)\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
		mismatch_set = true;
	}
	else
	{
		printf("\n  %sERROR%s: --mismatch [INT] has been set twice, please verify "
			"your choice\n\n", RED, COLOFF);
		printlist();
	}
} // ~Runopts::optMismatch

void Runopts::optGapOpen(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --gap_open [INT] requires a positive integer "
			"as input (ex. --gap_open 5)\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// set gap open
	if (!gap_open_set)
	{
		gap_open = atoi(argv[narg + 1]);
		if (gap_open < 0)
		{
			fprintf(stderr, "\n  %sERROR%s: --gap_open [INT] requires a positive "
				"integer as input (ex. --gap_open 5)\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
		gap_open_set = true;
	}
	else
	{
		printf("\n  %sERROR%s: --gap_open [INT] has been set twice, please verify "
			"your choice\n\n", RED, COLOFF);
		printlist();
	}
} // ~Runopts::optGapOpen

void Runopts::optGapExt(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --gap_ext [INT] requires a positive integer "
			"as input (ex. --gap_ext 2)\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// set gap extend
	if (!gap_ext_set)
	{
		gap_extension = atoi(argv[narg + 1]);
		if (gap_extension < 0)
		{
			fprintf(stderr, "\n  %sERROR%s: --gap_ext [INT] requires a positive "
				"integer as input (ex. --gap_ext 2)\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
		gap_ext_set = true;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: --gap_ext [INT] has been set twice, please "
			"verify your choice\n\n", RED, COLOFF);
		printlist();
	}
} // ~Runopts::optGapExt

void Runopts::optNumSeeds(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --num_seeds [INT] requires a positive integer "
			"as input (ex. --num_seeds 6)\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// set number of seeds
	if (seed_hits < 0)
	{
		char* end = 0;
		seed_hits = (int)strtol(argv[narg + 1], &end, 10); // convert to integer
		if (seed_hits <= 0)
		{
			fprintf(stderr, "\n  %sERROR%s: --num_seeds [INT] requires a positive "
				"integer (>0) as input (ex. --num_seeds 6)\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: --num_seeds [INT] has been set twice, please "
			"verify your choice\n\n", RED, COLOFF);
		printlist();
	}
} // ~Runopts::optNumSeeds

  /* --fastx */
void Runopts::optFastx(char **argv, int &narg)
{
	if (fastxout)
	{
		fprintf(stderr, "\n  %sERROR%s: --fastx has already been set once.\n\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		fastxout = true;
		narg++;
	}
} // ~Runopts::optFastx

void Runopts::optSam(char **argv, int &narg)
{
	if (samout)
	{
		fprintf(stderr, "\n  %sERROR%s: --sam has already been set once.\n\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		samout = true;
		narg++;
	}
} // ~Runopts::optSam

void Runopts::optBlast(char **argv, int &narg)
{
	std::stringstream ss;

	if (blastout)
	{
		ss << std::endl << "  " << RED << "ERROR" << COLOFF
			<< ": --blast [STRING] has already been set once."
			<< std::endl << std::endl;
		std::cerr << ss.str();
		exit(EXIT_FAILURE);
	}

	std::string str(argv[narg + 1]);
	// split blast options into vector by space
	std::istringstream iss(str);
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
			ss << std::endl << "  " << RED << "ERROR" << COLOFF
				<< ": `" << opt << "` is not supported in --blast [STRING]."
				<< std::endl << std::endl;
			std::cerr << ss.str();
			exit(EXIT_FAILURE);
		}
	}
	// more than 1 field with blast human-readable format given
	if (blast_human_readable && (blastops.size() > 1))
	{
		ss << std::endl << "  " << RED << "ERROR" << COLOFF
			<< ": for human-readable format, --blast [STRING] can only contain a single field '0'."
			<< std::endl << std::endl;
		std::cerr << ss.str();
		exit(EXIT_FAILURE);
	}
	// both human-readable and tabular format options have been chosen
	if (blast_human_readable && blastFormat == BlastFormat::TABULAR)
	{
		ss << std::endl << "  " << RED << "ERROR" << COLOFF
			<< ": --blast [STRING] can only have one of the options '0' (human-readable) or '1' (tabular)."
			<< std::endl << std::endl;
		std::cerr << ss.str();
		exit(EXIT_FAILURE);
	}

	blastout = true;
	narg += 2;
} // ~Runopts::optBlast

void Runopts::optMinLis(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --min_lis [INT] requires an integer (>=0) as "
			"input (ex. --min_lis 2) (note: 0 signifies to search all high scoring "
			"reference sequences).\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// min_lis_gv has already been set
	else if (min_lis_set)
	{
		fprintf(stderr, "\n  %sERROR%s: --min_lis [INT] has been set twice, please "
			"verify your choice.\n\n", RED, COLOFF);
		printlist();
	}
	else
	{
		if ((sscanf(argv[narg + 1], "%d", &min_lis) != 1) || (min_lis < 0))
		{
			fprintf(stderr, "\n  %sERROR%s: --min_lis [INT] must be >= 0 (0 signifies "
				"to search all high scoring reference sequences).\n\n",
				RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
		min_lis_set = true;
	}
} // ~Runopts::optMinLis

void Runopts::optBest(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --best [INT] requires an integer (> 0) "
			"as input (ex. --best 2).\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// best_set has already been set
	else if (best_set)
	{
		fprintf(stderr, "\n  %sERROR%s: --best [INT] has been set twice, please "
			"verify your choice.\n\n", RED, COLOFF);
		printlist();
	}
	else
	{
		if ((sscanf(argv[narg + 1], "%d", &num_best_hits) != 1))
		{
			fprintf(stderr, "\n  %sERROR%s: could not read --best [INT] as integer.\n\n",
				RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
		best_set = true;
	}
} // ~Runopts::optBest

  /* --num_alignments [INT] */
void Runopts::optNumAlignments(char **argv, int &narg)
{
	if (argv[narg + 1] == NULL)
	{
		fprintf(stderr, "\n  %sERROR%s: --num_alignments [INT] requires an integer "
			"(>=0) as input (ex. --num_alignments 2) (note: 0 signifies to output all alignments).\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// --num_alignments has already been set
	else if (num_alignments_set)
	{
		fprintf(stderr, "\n  %sERROR%s:--num_alignments [INT] has been set twice, please verify your command parameters.\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// set number of alignments to output reaching the E-value
	else
	{
		num_alignments = atoi(argv[narg + 1]);
		if (num_alignments < 0)
		{
			fprintf(stderr, "\n  %sERROR%s: --num_alignments [INT] must be >= 0 (0 signifies to output all alignments).\n\n",
				RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
		num_alignments_set = true;
	}
} // ~Runopts::optNumAlignments

void Runopts::optEdges(char **argv, int &narg)
{
	// --edges is already set
	if (edges_set)
	{
		fprintf(stderr, "\n  %sERROR%s: --edges [INT]%% has already been set once.\n\n",
			RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		char *end = 0;
		// find if % sign exists
		char* test = strpbrk(argv[narg + 1], "%");
		if (test != NULL)
			as_percent = true;
		// convert to integer
		edges = (int)strtol(argv[narg + 1], &end, 10);

		if (edges < 1 || edges > 10)
		{
			fprintf(stderr, "\n  %sERROR%s: --edges [INT]%% requires a positive integer "
				"between 0-10 as input (ex. --edges 4).\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}

		narg += 2;
	}
} // ~Runopts::optEdges

void Runopts::optFullSearch(char **argv, int &narg)
{
	if (full_search_set)
	{
		fprintf(stderr, "\n  %sERROR%s: BOOL --full_search has been set twice, please "
			"verify your choice.\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		full_search_set = true;
		full_search = true;
		narg++;
	}
} // ~Runopts::optFullSearch

void Runopts::optSQ(char **argv, int &narg)
{
	if (yes_SQ)
	{
		fprintf(stderr, "\n  %sERROR%s: BOOL --SQ has been set twice, please verify "
			"your choice.\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	else
	{
		yes_SQ = true;
		narg++;
	}
} // ~Runopts::optSQ

void Runopts::optPasses(char **argv, int &narg)
{
	if (passes_set)
	{
		fprintf(stderr, "\n  %sERROR%s: --passes [INT,INT,INT] has been set twice, "
			"please verify your choice.\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
	// set passes
	else
	{
		std::vector<uint32_t> skiplengths_v;
		char *end = 0;
		int32_t t = (int)strtol(strtok(argv[narg + 1], ","), &end, 10);
		if (t > 0) skiplengths_v.push_back(t);
		else
		{
			fprintf(stderr, "\n  %sERROR%s: all three integers in --passes [INT,INT,INT] "
				"must contain positive integers where 0<INT<(shortest read length)."
				"\n\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		t = (int)strtol(strtok(NULL, ","), &end, 10);
		if (t > 0) skiplengths_v.push_back(t);
		else
		{
			fprintf(stderr, "\n  %sERROR%s: all three integers in --passes [INT,INT,INT] "
				"must contain positive integers where 0<INT<(shortest read length). "
				"\n\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		t = (int)strtol(strtok(NULL, ","), &end, 10);
		if (t > 0) skiplengths_v.push_back(t);
		else
		{
			fprintf(stderr, "\n  %sERROR%s: all three integers in --passes [INT,INT,INT] "
				"must contain positive integers where 0<INT<(shortest read length)."
				"\n\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}

		skiplengths.push_back(skiplengths_v);
		narg += 2;
		passes_set = true;
	}
} // ~Runopts::optPasses

void Runopts::optId(char **argv, int &narg)
{
	// % id
	if (align_id < 0)
	{
		if ((sscanf(argv[narg + 1], "%lf", &align_id) != 1) ||
			(align_id < 0) || (align_id > 1))
		{
			fprintf(stderr, "\n  %sERROR%s: --id [DOUBLE] must be a positive float "
				"with value 0<=id<=1.\n\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: --id [DOUBLE] has been set twice, please "
			"verify your command parameters.\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
} // ~Runopts::optId

void Runopts::optCoverage(char **argv, int &narg)
{
	// % query coverage
	if (align_cov < 0)
	{
		if ((sscanf(argv[narg + 1], "%lf", &align_cov) != 1) ||
			(align_cov < 0) || (align_cov > 1))
		{
			fprintf(stderr, "\n  %sERROR%s: --coverage [DOUBLE] must be a positive "
				"float with value 0<=id<=1.\n\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: --coverage [DOUBLE] has been set twice, please "
			"verify your command parameters.\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
} // ~Runopts::optCoverage

void Runopts::optVersion(char **argv, int &narg)
{
	std::cout << std::endl
		<< "SortMeRNA version " << SORTMERNA_MAJOR << "." << SORTMERNA_MINOR << "." << SORTMERNA_PATCH << std::endl
		<< "Build Date: " << sortmerna_build_compile_date << std::endl
		<< sortmerna_build_git_sha << std::endl
		<< sortmerna_build_git_date << std::endl;
	exit(EXIT_SUCCESS);
} // ~Runopts::optVersion

void Runopts::optUnknown(char **argv, int &narg, char * opt)
{
	std::stringstream ss;
	ss << "\n  " << RED << "ERROR" << COLOFF << ": option --" << opt << " not recognized" << std::endl << std::endl;
	std::cout << ss.str();
	printlist();
} // ~Runopts::optUnknown

void Runopts::opt_e_Evalue(char **argv, int &narg)
{
	std::stringstream ss;

	// E-value
	if (argv[narg + 1] == NULL)
	{
		ss << "\n  " << RED << "ERROR" << COLOFF
			<< ": -e [DOUBLE] requires a positive double as input (ex. --e 1e-5)" << std::endl;
		std::cerr << ss.str();
		exit(EXIT_FAILURE);
	}

	if (evalue < 0)
	{
		sscanf(argv[narg + 1], "%lf", &evalue);
		if (evalue < 0)
		{
			fprintf(stderr, "\n  %sERROR%s: -e [DOUBLE] requires a positive double "
				"as input (ex. --e 1e-5)\n", RED, COLOFF);
			exit(EXIT_FAILURE);
		}
		narg += 2;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: -e [DOUBLE] has been set twice, please verify "
			"your command parameters.\n\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_e_Evalue

void Runopts::opt_F_ForwardOnly(char **argv, int &narg)
{
	// only forward strand
	if (!forward)
	{
		forward = true;
		narg++;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: BOOL -F has been set more than once, please check "
			"your command parameters.\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_F_ForwardOnly

void Runopts::opt_R_ReverseOnly(char **argv, int &narg)
{
	// only reverse strand
	if (!reverse)
	{
		reverse = true;
		narg++;
	}
	else
	{
		fprintf(stderr, "\n  %sERROR%s: BOOL -R has been set more than once, please check "
			"your command parameters.\n", RED, COLOFF);
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_R_ReverseOnly

void Runopts::opt_h_Help()
{
	welcome();
	printlist();
} // ~Runopts::opt_h_Help

void Runopts::opt_v_Verbose(int & narg)
{
	verbose = true;
	narg++;
} // ~Runopts::opt_v_Verbose

void Runopts::opt_N_MatchAmbiguous(char **argv, int &narg)
{
	// match ambiguous N's
	if (!match_ambiguous_N)
	{
		match_ambiguous_N = true;
		score_N = atoi(argv[narg + 1]);
		narg += 2;
	}
	else
	{
		std::stringstream ss;
		ss << std::endl << " " << RED << "ERROR" << COLOFF
			<< ": BOOL -N has been set more than once, please check your command parameters." << std::endl;
		std::cerr << ss.str();
		exit(EXIT_FAILURE);
	}
} // ~Runopts::opt_N_MatchAmbiguous

  /* Number Processor threads to use */
void Runopts::opt_a_numProcThreads(char **argv, int &narg)
{
	std::stringstream ss;

	if (argv[narg + 1] == NULL)
	{
		ss << "\n  " << RED << "ERROR" << COLOFF
			<< ": -a [INT] requires an integer for number of Processor threads (ex. -a 8)" << std::endl;
		std::cerr << ss.str();
		exit(EXIT_FAILURE);
	}

	num_proc_thread = atoi(argv[narg + 1]);
	narg += 2;
} // ~Runopts::opt_a_numProcThreads

  /* Number of threads to use */
void Runopts::opt_threads(char **argv, int &narg)
{
	std::stringstream ss;

	if (argv[narg + 1] == NULL)
	{
		ss << "\n  " << RED << "ERROR" << COLOFF
			<< ": --threads [INT:INT:INT] requires 3 integers for number of "
			<< "Read:Write:Processor threads (ex. --threads 1:1:8)" << std::endl;
		std::cerr << ss.str(); ss.str("");
		exit(EXIT_FAILURE);
	}

	std::istringstream strm(argv[narg + 1]);
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

	narg += 2;
} // ~Runopts::opt_threads


void Runopts::opt_threads_pp(char **argv, int &narg)
{
	std::stringstream ss;

	if (argv[narg + 1] == NULL)
	{
		ss << "\n  " << RED << "ERROR" << COLOFF
			<< ": --thpp [INT:INT] requires 2 integers for number of "
			<< "Read:Processor threads (ex. --thpp 1:1)" << std::endl;
		std::cerr << ss.str(); ss.str("");
		exit(EXIT_FAILURE);
	}

	std::istringstream strm(argv[narg + 1]);
	std::string tok;
	for (int i = 0; std::getline(strm, tok, ':'); ++i)
	{
		switch (i)
		{
		case 0: num_read_thread_pp = std::stoi(tok); break;
		case 1: num_proc_thread_pp = std::stoi(tok); break;
		}
	}

	narg += 2;
} // ~Runopts::opt_threads_pp


void Runopts::opt_threads_rep(char **argv, int &narg)
{
	std::stringstream ss;

	if (argv[narg + 1] == NULL)
	{
		ss << "\n  " << RED << "ERROR" << COLOFF
			<< ": --threp [INT:INT] requires 2 integers for number of "
			<< "Read:Processor threads (ex. --threp 1:1)" << std::endl;
		std::cerr << ss.str(); ss.str("");
		exit(EXIT_FAILURE);
	}

	std::istringstream strm(argv[narg + 1]);
	std::string tok;
	for (int i = 0; std::getline(strm, tok, ':'); ++i)
	{
		switch (i)
		{
		case 0: num_read_thread_rep = std::stoi(tok); break;
		case 1: num_proc_thread_rep = std::stoi(tok); break;
		}
	}

	narg += 2;
} // ~Runopts::opt_threads_rep

  // Required parameter
void Runopts::opt_d_KeyValDatabase(char **argv, int &narg)
{
	std::stringstream ss;

	// kvdb is specified
	if (argv[narg + 1] != NULL)
	{
		kvdbPath.assign(argv[narg + 1]); // e.g. "C:/a01_projects/clarity_genomics/data/kvdb"
		//ss << "\n  " << RED << "ERROR" << COLOFF
		//	<< ": -d [STRING] requires a folder path for key-value database (ex. -d /var/data/kvdb)" << std::endl;
		//std::cerr << ss.str();
		//exit(EXIT_FAILURE);
	}
	narg += 2;
} // ~Runopts::opt_d_KeyValDatabase


void Runopts::opt_debug_put_kvdb(int &narg)
{
	dbg_put_kvdb = true;
	narg++;
}


void Runopts::opt_Default(char **argv, int &narg)
{
	std::stringstream ss;
	ss << "\n  " << RED << "ERROR" << COLOFF << ": [Line " << __LINE__ << ": " << __FILE__
		<< "] '" << argv[narg][1] << "' is not one of the options." << std::endl;
	std::cerr << ss.str(); ss.str("");
	printlist();
} // ~Runopts::opt_Default

  /* Processing task */
void Runopts::optTask(char **argv, int &narg)
{
	int taskOpt = 4;
	sscanf(argv[narg + 1], "%d", &taskOpt);

	if (taskOpt > 4) {
		std::cerr << "Option −−task " << taskOpt << " Can only take values: [0..4] ... " << std::endl;
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

	narg += 2;
} // ~Runopts::optReport

  // interactive session '--cmd'
void Runopts::optInteractive(char **argv, int &narg)
{
	interactive = true;
	++narg;
} // ~Runopts::optInteractive

void Runopts::test_kvdb_path()
{
	if (kvdbPath.size() == 0)
	{
		kvdbPath = get_user_home() + "/kvdb";
		std::cout << __func__ << ": Key-value DB location was not specified. Setting default" << std::endl;
	}

	std::cout << __func__ << ": Using Key-value DB location: " << kvdbPath << std::endl;

	if (dirExists(kvdbPath))
	{
		// dir exists and not empty -> exception
		auto count = list_dir(kvdbPath);
		if (count > 0) 
		{
			if (ALIGN_REPORT::align == alirep || ALIGN_REPORT::all == alirep || ALIGN_REPORT::alipost == alirep)
			{
				std::cerr << __func__ << ": Directory " << kvdbPath
					<< " exists and is Not empty. Please, make sure the directory is empty, or specify a different directory using option '-d'" << std::endl;
				exit(1);
			}
		}
		else
		{
			if (ALIGN_REPORT::postproc == alirep || ALIGN_REPORT::report == alirep)
			{
				std::cerr << __func__ << ": Directory " << kvdbPath
					<< " is empty. Alignment has to be performed first. Please, use option '--task 0 | 3 | 4'" << std::endl;
				exit(1);
			}
			// dir exists and empty -> use
		}
	}
	else
	{
		// dir does not exist -> try creating
		std::cout << __func__ << ": Directory " << kvdbPath << " does not exists - will attempt to create";
	}
} // ~test_kvdb_path


void Runopts::process(int argc, char**argv, bool dryrun)
{
	if (dryrun) return;

	int narg = 1; // parse the command line input

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
	//else
	//{
	//	maxpages_gv = size / pagesize_gv;
	//}
#else
				  //maxpages_gv = sysconf(_SC_PHYS_PAGES);
#endif

#if defined(_WIN32)
	_setmode(_fileno(stderr), _O_BINARY);
#endif

	if (argc == 1)
	{
		verbose = true;
		welcome();
		fprintf(stderr, "  For help or more information on usage, type `./sortmerna %s-h%s'\n\n", BOLD, COLOFF);
		exit(EXIT_FAILURE);
	}

	while (narg < argc)
	{
		switch (argv[narg][1])
		{
			// options beginning with '--'
		case '-':
		{
			char* opt = argv[narg];
			opt += 2; // skip the '--'

					  // FASTA/FASTQ reads sequences
			if (strcmp(opt, "reads") == 0) optReads(argv, narg);
#ifdef HAVE_LIBZ
			// FASTA/FASTQ compressed reads sequences
			else if (strcmp(opt, "reads-gz") == 0) optReadsGz(argv, narg);
#endif
			// FASTA reference sequences
			else if (strcmp(opt, "ref") == 0) optRef(argv, narg);
			// the name of output aligned reads
			else if (strcmp(opt, "aligned") == 0) optAligned(argv, narg);
			// the name of output rejected reads
			else if (strcmp(opt, "other") == 0) optOther(argv, narg);
			// output overall statistics file
			else if (strcmp(opt, "log") == 0) optLog(argv, narg);
			// output FASTA/FASTQ reads passing E-value threshold but having < %id 
			// and < %coverage scores for de novo OTU construction
			else if (strcmp(opt, "de_novo_otu") == 0) optDeNovoOtu(argv, narg);
			// output OTU map
			else if (strcmp(opt, "otu_map") == 0) optOtuMap(argv, narg);
			// output non-aligned reads to SAM/BLAST files
			else if (strcmp(opt, "print_all_reads") == 0) optPrintAllReads(argv, narg);
			// don't add pid to output files
			else if (strcmp(opt, "pid") == 0) optPid(argv, narg);
			// put both paired reads into --accept reads file
			else if (strcmp(opt, "paired_in") == 0) optPairedIn(argv, narg);
			// put both paired reads into --other reads file
			else if (strcmp(opt, "paired_out") == 0) optPairedOut(argv, narg);
			// the score for a match
			else if (strcmp(opt, "match") == 0) optMatch(argv, narg);
			// the score for a mismatch
			else if (strcmp(opt, "mismatch") == 0) optMismatch(argv, narg);
			// the score for a gap
			else if (strcmp(opt, "gap_open") == 0) optGapOpen(argv, narg);
			// the score for a gap extension
			else if (strcmp(opt, "gap_ext") == 0) optGapExt(argv, narg);
			// number of seed hits before searching for candidate LCS
			else if (strcmp(opt, "num_seeds") == 0) optNumSeeds(argv, narg);
			// output all hits in FASTX format
			else if (strcmp(opt, "fastx") == 0) optFastx(argv, narg);
			// output all hits in SAM format
			else if (strcmp(opt, "sam") == 0) optSam(argv, narg);
			// output all hits in BLAST format
			else if (strcmp(opt, "blast") == 0) optBlast(argv, narg);
			// output best alignment as predicted by the longest increasing subsequence
			else if (strcmp(opt, "min_lis") == 0) optMinLis(argv, narg);
			// output best alignment as predicted by the longest increasing subsequence
			else if (strcmp(opt, "best") == 0) optBest(argv, narg);
			// output all alignments
			else if (strcmp(opt, "num_alignments") == 0) optNumAlignments(argv, narg);
			// number of nucleotides to add to each edge of an alignment region before extension
			else if (strcmp(opt, "edges") == 0) optEdges(argv, narg);
			// execute full index search for 0-error and 1-error seed matches
			else if (strcmp(opt, "full_search") == 0) optFullSearch(argv, narg);
			// do not output SQ tags in the SAM file
			else if (strcmp(opt, "SQ") == 0) optSQ(argv, narg);
			else if (strcmp(opt, "passes") == 0) optPasses(argv, narg); // --passes
			else if (strcmp(opt, "id") == 0) optId(argv, narg);
			else if (strcmp(opt, "coverage") == 0) optCoverage(argv, narg);
			else if (strcmp(opt, "version") == 0) optVersion(argv, narg); // version number
			else if (strcmp(opt, "task") == 0) optTask(argv, narg);
			else if (strcmp(opt, "cmd") == 0) optInteractive(argv, narg); // '--cmd' interactive session
																		  // threads
			else if (strcmp(opt, "thread") == 0) opt_threads(argv, narg); // '--thread 1:1:8' num alignment threads
			else if (strcmp(opt, "thpp") == 0) opt_threads_pp(argv, narg); // '--thpp 1:1' num post-proc threads
			else if (strcmp(opt, "threp") == 0) opt_threads_rep(argv, narg); // '--threp 1:1' num report threads
			else if (strcmp(opt, "dbg_put_db") == 0) opt_debug_put_kvdb(narg); // '--dbg_put_db'
			else optUnknown(argv, narg, opt);
		}
		break;
		case 'a': opt_a_numProcThreads(argv, narg); break;
		case 'd': opt_d_KeyValDatabase(argv, narg); break; // required
		case 'e': opt_e_Evalue(argv, narg);	break;
		case 'F': opt_F_ForwardOnly(argv, narg); break;
		case 'R': opt_R_ReverseOnly(argv, narg); break;
		case 'h': opt_h_Help(); break;
		case 'v': opt_v_Verbose(narg); break;
		case 'N': opt_N_MatchAmbiguous(argv, narg); break;
		default: opt_Default(argv, narg);
		}//~switch
	}//~while ( narg < argc )

	// validate the options
	test_kvdb_path();

	 // ERROR messages ******* 
	 // Reads file is mandatory
	if (readsfile.empty() || indexfiles.empty())
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] a reads file (--reads file.{fa/fq}) and a "
			"reference sequence file (--ref /path/to/file1.fasta,/path/to/index1) "
			"are mandatory input.\n\n", RED, COLOFF, __LINE__, __FILE__);
		printlist();
	}

	// Basename for aligned reads is mandatory
	if (filetype_ar.size() == 0)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] parameter --aligned [STRING] is mandatory.\n\n",
			RED, COLOFF, __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	// No output format has been chosen
	else if (!(fastxout || blastout || samout || otumapout || doLog || de_novo_otu))
	{
		fprintf(stderr,
			"\n  %sERROR%s: [Line %d: %s] no output format has been chosen (fastx/sam/blast/otu_map/log).\n\n",
			RED, COLOFF, __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	// Options --paired_in and --paired_out can only be used with FASTA/Q output
	if (!fastxout && (pairedin || pairedout))
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] options --paired_in and --paired_out "
			"must be accompanied by option --fastx.\n", RED, COLOFF, __LINE__, __FILE__);
		fprintf(stderr, "  These BOOLs are for FASTA and FASTQ output files, for "
			"maintaining paired reads together.\n");
		exit(EXIT_FAILURE);
	}

	// Basename for non-aligned reads is mandatory
	if (filetype_or.size() != 0)
	{
		if (!fastxout && (blastout || samout))
		{
			fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] option --other [STRING] can only be used together "
				"with the --fastx option.\n\n", RED, COLOFF, __LINE__, __FILE__);
			exit(EXIT_FAILURE);
		}
	}

	// An OTU map can only be constructed with the single best alignment per read
	if (otumapout && num_alignments_set)
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --otu_map cannot be set together with "
			"--num_alignments [INT].\n", RED, COLOFF, __LINE__, __FILE__);
		fprintf(stderr, "  The option --num_alignments [INT] doesn't keep track of "
			"the best alignment which is required for constructing an OTU map.\n");
		fprintf(stderr, "  Use --otu_map with --best [INT] instead.\n\n");
		exit(EXIT_FAILURE);
	}

	// If --num_alignments output was chosen, check an alignment format has also been chosen
	if (num_alignments_set && !(blastout || samout || fastxout))
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --num_alignments [INT] has been set but no output "
			"format has been chosen (--blast, --sam or --fastx).\n\n", RED, COLOFF, __LINE__, __FILE__);
		exit(EXIT_FAILURE);
	}

	// If --best output was chosen, check an alignment format has also been chosen
	if (best_set && !(blastout || samout || otumapout))
	{
		fprintf(stderr, "\n  %sERROR%s: [Line %d: %s] --best [INT] has been set but no output "
			"format has been chosen (--blast or --sam or --otu_map).\n\n", RED, COLOFF, __LINE__, __FILE__);
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
	if (verbose) welcome();
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
		if (fastxout && !(blastout || samout || otumapout || doLog || de_novo_otu))
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

	for (int i = 0; i < argc; i++) {
		cmdline.append(argv[i]);
		cmdline.append(" ");
	}
} // ~Runopts::process




  /*! @fn welcome()
  *  @brief outputs copyright, disclaimer and contact information
  *  @param none
  #  @return none
  */
void welcome()
{
	std::stringstream ss;

	ss << std::endl
		<< "  Program:      SortMeRNA version " << SORTMERNA_MAJOR << "." << SORTMERNA_MINOR << "." << SORTMERNA_PATCH << std::endl
		<< "  Copyright:    2016-2018 Clarity Genomics BVBA:" << std::endl
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
void printlist()
{
	std::stringstream ss;

	ss << std::endl
		<< "  usage:   ./sortmerna --ref db.fasta,db.idx --reads file.fa --aligned base_name_output [OPTIONS]:" << std::endl
#ifdef HAVE_LIBZ
		<< "  OR" << std::endl
		<< "  usage:   ./sortmerna --ref db.fasta,db.idx --reads-gz file.fa.gz --aligned base_name_output [OPTIONS]:" << std::endl
#endif
		<< std::endl
		<< "  -------------------------------------------------------------------------------------------------------------"  << std::endl
		<< "  | option              type-format       description                                              default    |"  << std::endl
		<< "  -------------------------------------------------------------------------------------------------------------"  << std::endl
		<< "  [REQUIRED OPTIONS]: "                                                                                           << std::endl << BOLD
		<< "    --ref            "                                                                                            << COLOFF << UNDL
		<<                       "  STRING,STRING"                                                                            << COLOFF
		<<                                       "   FASTA reference file:index file                           "              << GREEN 
		<<                                                                                                     "mandatory"    << COLOFF << std::endl
		<< "                                         If passing multiple reference files, separate them"                      << std::endl
		<< "                                         using the delimiter ':' (Linux) or ';' (Windows),"                       << std::endl
		<< "                 (ex. --ref /path/to/file1.fasta,/path/to/index1:/path/to/file2.fasta,path/to/index2)"            << std::endl << BOLD
		<< "    --reads          "                                                                                            << COLOFF << UNDL
		<<                       "  STRING       "                                                                            << COLOFF
		<<                                       "   FASTA/FASTQ raw reads file                                "              << GREEN 
		<<                                                                                                     "mandatory"    << COLOFF << std::endl
#ifdef HAVE_LIBZ
		<< "        OR"                                                                                                       << std::endl << BOLD                                                                                                    
		<< "    --reads-gz       "                                                                                            << COLOFF << UNDL
		<<                       "  STRING       "                                                                            << COLOFF
		<<                                       "   FASTA/FASTQ compressed (with gzip) reads file             "              << GREEN 
		<<                                                                                                     "mandatory"    << COLOFF << std::endl << BOLD
#endif
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

	exit(EXIT_FAILURE);
}//~printlist()