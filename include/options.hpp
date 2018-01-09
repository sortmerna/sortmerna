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
	std::string kvdbPath = "C:/a01_projects/clarity_genomics/data/kvdb"; // key-value database for alignment results
	std::string readsfile; // '--reads | --reads-gz' reads file path
	std::string filetype_ar; // '--aligned' aligned reads output file
	std::string filetype_or; // '--other' rejected reads output file
	std::string cmdline;

	int numcpu = -1; // '-a' number of threads to use TODO: remove (see num_proc_threads)
	int num_fread_threads = 1; // number of threads reading the Reads file.
	int num_proc_threads = 4;  // '-a' number of threads to use for processing

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
	int32_t min_lis = -1; // '--min_lis' search all alignments having the first INT longest LIS
	int32_t seed_hits = -1;
	int32_t num_best_hits = 0;
	int32_t edges = -1;

	bool forward = false; // '-F' search only the forward strand
	bool reverse = false; // '-R' search only the reverse-complementary strand
	bool pairedin = false; // '--paired_in' both paired-end reads go in --aligned fasta/q file (interleaved reads only, see Section 4.2.4 of User Manual)
	bool pairedout = false; // '--paired_out' both paired-end reads go in --other fasta/q file (interleaved reads only, see Section 4.2.4 of User Manual)
	bool de_novo_otu = false; // '--de_novo_otu' FASTA/FASTQ file for reads matching database < %%id (set using --id) and < %%cov (set using --coverage)
	bool doLog = false; // '--log' output overall statistics
	bool print_all_reads = false; // '--print_all_reads' output null alignment strings for non-aligned reads to SAM and/or BLAST tabular files
	bool samout = false; // '--sam' output SAM alignment (for aligned reads only)
	bool blastout = false; // '--blast' output alignments in various Blast-like formats
	bool fastxout = false; // '--fastx' output FASTA/FASTQ file (for aligned and/or rejected reads)
	bool otumapout = false; // '--otu_map' output OTU map (input to QIIME's make_otu_table.py)
	bool map_size_set = false; // TODO: remove - mmap related - not used
	bool pid = false;
	bool as_percent = false;
	bool full_search = false;
	bool exit_early = false; // flag to exit processing when either the reads or the reference file is empty or not FASTA/FASTQ
	bool verbose; /* Verbose mode */
	bool have_reads_gz = false; // '--reads-gz' flags reads file is compressed and can be read
	bool yes_SQ = false; // --SQ add SQ tags to the SAM file

	enum ALIGN_REPORT { align, postproc, report, all };
	ALIGN_REPORT alirep = align;
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
	/* '--refs' (required): Pairs (Reference file, Index name) */
	std::vector<std::pair<std::string, std::string>> indexfiles;
	/* '--passes' (optional): for each index file three intervals at which to place the seed on the read. <-- Refstats::load
		Defaults: 0 */
	std::vector<std::vector<uint32_t>> skiplengths;

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
	bool match_ambiguous_N = false;
	bool min_lis_set = false;
	bool num_alignments_set = false;
	bool best_set = false;
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
	void opt_e_Evalue(char **argv, int &narg);
	void opt_a_NumCpus(char **argv, int &narg);
	void opt_F_ForwardOnly(char **argv, int &narg);
	void opt_R_ReverseOnly(char **argv, int &narg);
	void opt_h_Help();
	void opt_v_Verbose(int & narg);
	void opt_N_MatchAmbiguous(char **argv, int &narg);
	void opt_Default(char **argv, int &narg);
};
