#ifndef GUMBEL_PARAMS___HPP
#define GUMBEL_PARAMS___HPP

/* $Id: $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: gumbel_params.hpp

Author: Sergey Sheetlin

Contents: Wrapper classes for real time Gumbel parameters computing code

******************************************************************************/

#include <vector>
using namespace std;

#include "sls_pvalues.hpp"

class CGumbelParamsCalc 
{
public:

	static void Run2(
		long int rand_,//randomization number
		string randout_,//if defined, then the program outputs complete randomization information into a file

		long int gapopen_,//gap opening penalty
		long int gapopen1_,//gap opening penalty for a gap in the sequence #1
		long int gapopen2_,//gap opening penalty for a gap in the sequence #2

		long int gapextend_,//gap extension penalty
		long int gapextend1_,//gap extension penalty for a gap in the sequence #1
		long int gapextend2_,//gap extension penalty for a gap in the sequence #2

		string scoremat_file_name_,//scoring matrix file name
		string freqs1_file_name_,//probabilities1 file name
		string freqs2_file_name_,//probabilities1 file name
		double max_time_,//maximum allowed calculation time in seconds
		double max_mem_,//maximum allowed memory usage in MB
		double eps_lambda_,//relative error for lambda calculation
		double eps_K_,//relative error for K calculation
		string gumbelparout_file_name_,
		bool gapped_,

		string gumbelparin_file_name_,//Gumbel parameters input file name

		long int seqlen1_,//length of the sequence 1
		long int seqlen2_,//length of the sequence 2

		long int score1_,
		long int score2_,//P-values are calculated in the range [score1_,score2_]

		string pvalout_file_name_,//P-values file name
		bool insertions_after_deletions_,//if true, then insertions after deletions are allowed

		bool Gumbel_mode_,//true if Gumbel parameters will be calculated
		bool Pvalues_mode_);//true if P-values will be calculated




	//calculation of Gumbel parameters
	static void Params_Run2(
		long int rand_,//randomization number
		string randout_,//if defined, then the program outputs complete randomization information into a file

		long int open_,//gap opening penalty
		long int open1_,//gap opening penalty for a gap in the sequence #1
		long int open2_,//gap opening penalty for a gap in the sequence #2


		long int epen_,//gap extension penalty
		long int epen1_,//gap extension penalty for a gap in the sequence #1
		long int epen2_,//gap extension penalty for a gap in the sequence #2

		string smatr_file_name_,//scoring matrix file name
		string RR1_file_name_,//probabilities1 file name
		string RR2_file_name_,//probabilities1 file name
		double max_time_,//maximum allowed calculation time in seconds
		double max_mem_,//maximum allowed memory usage in MB
		double eps_lambda_,//relative error for lambda calculation
		double eps_K_,//relative error for K calculation
		string out_file_name_,
		bool gapped_,
		bool insertions_after_deletions_,//if true, then insertions after deletions are allowed
		Sls::set_of_parameters &gumbel_params_,
		long int match,//NEW - SW score for a match,
		long int mismatch,//NEW - SW score for a mismatch,
    double A_,//NEW - background frequency for A
		double C_,//NEW - background frequency for C
		double G_,//NEW - background frequency for G
		double T_,//NEW- background frequency for T
		double &lambda,
		double &K
		);

	static void Output_Params(
		Sls::set_of_parameters &gumbel_params_,
		string gumbelparout_file_name_);

	static void Read_Params(
		Sls::set_of_parameters &gumbel_params_,
		string gumbelparin_file_name_);



	//calculation of P-values

	static void P_values_Run2(
			long int Score1_,
			long int Score2_,
			long int Seq1Len_,
			long int Seq2Len_,
			Sls::set_of_parameters &pars_,
			vector<double>& pv_,
			vector<double>& pv_err_);

	static void Output_Pvalues(
		long int score1_,
		long int score2_,
		vector<double> &pv_,
		vector<double> &pv_err_,
		string pvalout_file_name_);//P-values file name




};

#endif // GUMBEL_PARAMS___HPP


