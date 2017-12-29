#pragma once
/**
 * FILE: options.hpp
 * Created: Aug 19, 2017 Sat
 *
 * skiplength
 *    skip lengths for pass 1, pass 2 and pass 3 in first step of sortmerna
 *    pipeline for each reference database searched
 */
#include <vector>
#include <string>

#include "common.hpp"

struct Runopts {
	// required '--refs': Pairs (Reference file, Index name)
	std::vector<std::pair<std::string, std::string>> indexfiles;
	// optional '--passes': for each index file three intervals at which to place the seed on the read.
	// init to 0 if '--passes' not provided. Set in Index::load_stats
	std::vector<std::vector<uint32_t>> skiplengths;
	std::string kvdbPath = "C:/a01_projects/clarity_genomics/data/kvdb"; // key-value database for match results

	int num_cpus = 1; // number of CPUs on this machine
	int num_fread_threads = 1; // number of threads reading the Reads file.
	int num_proc_threads = 4;  // '-a' number of threads to use for processing
	std::string readsfile; // '--reads | --reads-gz' reads file path
	char* ptr_filetype_ar = 0; // '--aligned' aligned reads output file
	char* ptr_filetype_or = 0; // '--other' rejected reads output file
	double evalue = -1.0; // '-e' E-value threshold
	long match = 2; // '--match' SW score (positive integer) for a match               TODO: change to int8_t
	long mismatch = -3; // '--mismatch' SW penalty (negative integer) for a mismatch   TODO: change to int8_t
	long gap_open = 5; // '--gap_open' SW penalty (positive integer) for introducing a gap
	long gap_extension = 2; // '--gap_ext' SW penalty (positive integer) for extending a gap
	long score_N = 0; // '-N' SW penalty for ambiguous letters (N's)                   TODO: change to int8_t
	bool exit_early = false; // flag to exit processing when either the reads or the reference file is empty or not FASTA/FASTQ
	double align_cov = -1.0; // '--coverage' query coverage threshold (the alignment must still pass the E-value threshold)

	// outputs the first '--num_alignments' found, unlike '--best', which searches
	// many alignments(specified by '--min_lis') prior to outputting the best ones.
	int32_t num_alignments = -1; // '--num_alignments'
	bool forward = false; // '-F' search only the forward strand
	bool reverse = false; // '-R' search only the reverse-complementary strand
	int numcpu = -1; // '-a' number of threads to use TODO: remove (see num_proc_threads)
	bool pairedin = false; // '--paired_in' both paired-end reads go in --aligned fasta/q file (interleaved reads only, see Section 4.2.4 of User Manual)
	bool pairedout = false; // '--paired_out' both paired-end reads go in --other fasta/q file (interleaved reads only, see Section 4.2.4 of User Manual)
	bool de_novo_otu = false; // '--de_novo_otu' FASTA/FASTQ file for reads matching database < %%id (set using --id) and < %%cov (set using --coverage)
	bool doLog = false; // '--log' output overall statistics
	bool print_all_reads = false; // '--print_all_reads' output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files
	bool samout = false; // '--sam' output SAM alignment (for aligned reads only)
	bool blastout = false; // '--blast' output alignments in various Blast-like formats
	bool fastxout = false; // '--fastx' output FASTA/FASTQ file (for aligned and/or rejected reads)
	bool otumapout = false; // '--otu_map' output OTU map (input to QIIME's make_otu_table.py)
	int32_t min_lis = -1; // '--min_lis' search all alignments having the first INT longest LIS
	BlastFormat blastFormat = BlastFormat::TABULAR;

	std::string cmdline;

	enum ALIGN_REPORT { align, report, both };
	ALIGN_REPORT alirep = align;

	bool have_reads_gz = false; // '--reads-gz' flags reads file is compressed and can be read
	bool yes_SQ = false; // --SQ add SQ tags to the SAM file

	Runopts(int argc, char**argv, bool dryrun)
	{ 
		process(argc, argv, dryrun);
		if (skiplengths.empty())
		{
			for ( int i = 0; i < indexfiles.size(); ++i )
			{
				skiplengths.push_back({ 0,0,0 });
			}
		}
	}
	~Runopts() {}

private:
	// SW alignment parameters
	bool match_set = false;
	bool mismatch_set = false;
	bool gap_open_set = false;
	bool gap_ext_set = false;
	bool full_search_set = false;
	bool passes_set = false;
	bool edges_set = false;
	bool match_ambiguous_N_gv = false;
	bool min_lis_gv_set = false;
	bool num_alignments_gv_set = false;
	bool best_gv_set = false;
	bool have_reads = false; // '--reads' flags reads file is plain text and can be read

private:
	// Functions
	void process(int argc, char**argv, bool dryrun);
	void optReads(char **argv, int &narg);
	void optReadsGz(char **argv, int &narg);
	void optRef(char **argv, int &narg);
	void optAligned(char **argv, int &narg);
	void optOther(char **argv, int &narg);
	void optLog(char **argv, int &narg);
	void optDeNovoOtu(char **argv, int &narg);
	void optOtuMap(char **argv, int &narg);
	void optPrintAllReads(char **argv, int &narg);
	void optPid(char **argv, int &narg);
	void optPairedIn(char **argv, int &narg);
	void optPairedOut(char **argv, int &narg);
	void optMatch(char **argv, int &narg);
	void optMismatch(char **argv, int &narg);
	void optGapOpen(char **argv, int &narg);
	void optGapExt(char **argv, int &narg);
	void optNumSeeds(char **argv, int &narg);
	void optFastx(char **argv, int &narg);
	void optSam(char **argv, int &narg);
	void optBlast(char **argv, int &narg);
	void optMinLis(char **argv, int &narg);
	void optBest(char **argv, int &narg);
	void optNumAlignments(char **argv, int &narg);
	void optEdges(char **argv, int &narg);
	void optFullSearch(char **argv, int &narg);
	void optSQ(char **argv, int &narg);
	void optPasses(char **argv, int &narg);
	void optId(char **argv, int &narg);
	void optCoverage(char **argv, int &narg);
	void optVersion(char **argv, int &narg);
	void optReport(char **argv, int &narg);
	void optUnknown(char **argv, int &narg, char * opt);
	void opt_a_NumCpus(char **argv, int &narg);
	void opt_e_Evalue(char **argv, int &narg);
	void opt_F_ForwardOnly(char **argv, int &narg);
	void opt_R_ReverseOnly(char **argv, int &narg);
	void opt_h_Help();
	void opt_v_Verbose(int & narg);
	void opt_N_MatchAmbiguous(char **argv, int &narg);
	void opt_m_MemMapSize(char **argv, int &narg);
	void opt_Default(char **argv, int &narg);
};
