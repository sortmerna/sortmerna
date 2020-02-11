/*
 * FILE: refstats.cpp
 * Created: Dec 23, 2017 Sat
 * @copyright 2016-20 Clarity Genomics BVBA
 */
#include <utility>
#include <fstream> // ifstream
#include <sstream>
#include <ios>
#include <vector>

#include "sls_alignment_evaluer.hpp" // ../alp/

#include "refstats.hpp"
#include "readstats.hpp"
#include "options.hpp"
#include "indexdb.hpp"

Refstats::Refstats(Runopts & opts, Readstats & readstats)
	:
	num_index_parts(opts.indexfiles.size(), 0),
	full_ref(opts.indexfiles.size(), 0),
	full_read(opts.indexfiles.size(), readstats.all_reads_len),
	lnwin(opts.indexfiles.size(), 0),
	partialwin(opts.indexfiles.size(), 0),
	minimal_score(opts.indexfiles.size(), 0),
	gumbel(opts.indexfiles.size(), std::pair<double, double>(-1.0, -1.0)),
	numbvs(opts.indexfiles.size(), 0),
	numseq(opts.indexfiles.size(), 0)
{
	std::stringstream ss;
	ss << STAMP << "Index Statistics calculation Start ...";
	std::cout << ss.str();

	auto starts = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed;

	load(opts, readstats);

	elapsed = std::chrono::high_resolution_clock::now() - starts;
	ss.str("");
	ss << STAMP << "Done. Time elapsed: " << std::setprecision(2) << std::fixed << elapsed.count() << " sec" << std::endl;
	std::cout << ss.str();
}

/**
 * load reference statistics stored in the '.stats' files 
 */
void Refstats::load(Runopts & opts, Readstats & readstats)
{
	std::stringstream ss;

	// create and initialize scoring matrix
	long alphabetSize = 4;
	long **scoring_matrix = new long *[alphabetSize];

	for (long i = 0; i < alphabetSize; i++)
	{
		scoring_matrix[i] = new long[alphabetSize];
		for (long j = 0; j < alphabetSize; j++)
		{
			if (i == j)
				scoring_matrix[i][j] = opts.match;
			else
				scoring_matrix[i][j] = opts.mismatch;
		}
	}

	// loop through the .stats index files for each database
	for (uint16_t index_num = 0; index_num < (uint16_t)opts.indexfiles.size(); index_num++)
	{
		std::ifstream stats(opts.indexfiles[index_num].second + ".stats", std::ios::in | std::ios::binary);
		if (!stats.good())
		{
			// should never come here. Index is built and validated prior this call.
			ss.str("");
			ss << STAMP << "Cannot open the index file [" << opts.indexfiles[index_num].second << ".stats]";
			ERR(ss.str());
			exit(EXIT_FAILURE);
		}

		// read the file size for file used to build the index
		size_t filesize = 0;
		stats.read(reinterpret_cast<char*>(&filesize), sizeof(size_t));

		// read the fasta file name used to build the index
		uint32_t fastafile_len = 0;
		stats.read(reinterpret_cast<char*>(&fastafile_len), sizeof(fastafile_len));
		char fastafile_name[2000];
		stats.read(reinterpret_cast<char*>(fastafile_name), sizeof(char)*fastafile_len);

		// compute reference database file size for this index
		//FILE *fastafile = ::fopen(opts.indexfiles[index_num].first.c_str(), "r");
		//if (fastafile == NULL)
		//{
		//	ss.str("");
		//	ss <<  "could not open FASTA reference file: " << opts.indexfiles[index_num].first;
		//	ERR(ss.str());
		//	exit(EXIT_FAILURE);
		//}

		// A,C,G,T background frequencies to compute the Gumbel parameters lambda and K
		double background_freq_gv[4] = { 0 };
		// A/C/G/T distribution frequencies
		stats.read(reinterpret_cast<char*>(&background_freq_gv), sizeof(double) * 4);
		// total length of sequences in the complete database
		stats.read(reinterpret_cast<char*>(&full_ref[index_num]), sizeof(uint64_t));
		// sliding window length lnwin & initialize
		stats.read(reinterpret_cast<char*>(&lnwin[index_num]), sizeof(uint32_t));
		// total number of reference sequences in one complete reference database
		stats.read(reinterpret_cast<char*>(&numseq[index_num]), sizeof(uint64_t));
		partialwin[index_num] = lnwin[index_num] / 2;
		// number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward
		numbvs[index_num] = 4 * (partialwin[index_num] - 3);

		// set the window shift for different seed lengths (if not set by user, or one of the lengths is <= 0)
		if ((opts.skiplengths[index_num][0] == 0) || (opts.skiplengths[index_num][1] == 0) || (opts.skiplengths[index_num][2] == 0))
		{
			opts.skiplengths[index_num][0] = lnwin[index_num];
			opts.skiplengths[index_num][1] = partialwin[index_num];
			opts.skiplengths[index_num][2] = 3;
		}

		// get number of index parts i.e. how many parts the index has
		stats.read(reinterpret_cast<char*>(&num_index_parts[index_num]), sizeof(uint16_t));

		std::vector<index_parts_stats> hold;

		// information on the location and size of sequences used to build each index part
		for (uint16_t j = 0; j < num_index_parts[index_num]; j++)
		{
			index_parts_stats stats_hold;
			stats.read(reinterpret_cast<char*>(&stats_hold), sizeof(index_parts_stats));
			hold.push_back(stats_hold);
		}

		index_parts_stats_vec.push_back(hold);

		// Gumbel parameters
		long **substitutionScoreMatrix = scoring_matrix;
		long gapOpen1 = opts.gap_open;
		long gapOpen2 = opts.gap_open;
		long gapEpen1 = opts.gap_extension;
		long gapEpen2 = opts.gap_extension;
		bool insertions_after_deletions = false;
		double max_time = -1; // required if randomization parameters are set
		double max_mem = 500;
		double eps_lambda = 0.001;
		double eps_K = 0.005;
		long randomSeed = 182345345;
		double *letterFreqs1 = new double[alphabetSize];
		double *letterFreqs2 = new double[alphabetSize];
		long number_of_samples = 14112;
		long number_of_samples_for_preliminary_stages = 39;

		for (long i = 0; i < alphabetSize; i++)
		{
			// background probabilities for ACGT based on reference file
			letterFreqs1[i] = background_freq_gv[i];
			letterFreqs2[i] = background_freq_gv[i];
		}

		Sls::AlignmentEvaluer gumbelCalculator; // object to store the Gumbel parameters

		// set the randomization parameters
		// (will yield the same Lamba and K values on subsequent runs with the same input files)
		gumbelCalculator.set_gapped_computation_parameters_simplified(
			max_time,
			number_of_samples,
			number_of_samples_for_preliminary_stages);

		gumbelCalculator.initGapped(
			alphabetSize,
			substitutionScoreMatrix,
			letterFreqs1,
			letterFreqs2,
			gapOpen1,
			gapEpen1,
			gapOpen2,
			gapEpen2,
			insertions_after_deletions,
			eps_lambda,
			eps_K,
			max_time,
			max_mem,
			randomSeed);

		gumbel[index_num].first = gumbelCalculator.parameters().lambda;
		gumbel[index_num].second = gumbelCalculator.parameters().K;

		delete[] letterFreqs2;
		delete[] letterFreqs1;

		// Shannon's entropy for reference sequence nucleotide distribution
		double entropy_H_gv =
			-(background_freq_gv[0] * (log(background_freq_gv[0]) / log(2))
				+ background_freq_gv[1] * (log(background_freq_gv[1]) / log(2))
				+ background_freq_gv[2] * (log(background_freq_gv[2]) / log(2))
				+ background_freq_gv[3] * (log(background_freq_gv[3]) / log(2)));

		// Length correction for Smith-Waterman alignment score
		uint64_t expect_L = static_cast<uint64_t>(log((gumbel[index_num].second)*full_read[index_num] * full_ref[index_num]) / entropy_H_gv);

		// correct the reads & databases sizes for E-value calculation
		if (full_ref[index_num] > (expect_L*numseq[index_num]))
			full_ref[index_num] -= (expect_L*numseq[index_num]);

		full_read[index_num] -= (expect_L * readstats.all_reads_count);

		// minimum score required to reach E-value 
		// S = ln(E/Kmn)/-λ   <--   E = K*m*n*exp(-λS)
		minimal_score[index_num] = static_cast<uint32_t>(
			(log(opts.evalue
				/ ((double)(gumbel[index_num].second)
					* full_ref[index_num]
					* full_read[index_num])))
			/ -(gumbel[index_num].first));

		stats.close();
	} // ~for loop indices

	// free memory
	for (long i = 0; i < alphabetSize; i++)
	{
		delete[] scoring_matrix[i];
	};

	delete[] scoring_matrix;
} // ~Index::load_stats