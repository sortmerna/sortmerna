#ifndef INCLUDED_SLS_PVALUES
#define INCLUDED_SLS_PVALUES

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

File name: sls_pvalues.hpp

Author: Sergey Sheetlin

Contents: P-values calculation routines

******************************************************************************/


#include <vector>
#include <string>

#include "sls_normal_distr_array.hpp"


namespace Sls {

	struct set_of_parameters
	{
		double lambda;
		double lambda_error;

		double C;
		double C_error;


		double K;
		double K_error;

		double a_I;
		double a_I_error;

		double a_J;
		double a_J_error;

		double sigma;
		double sigma_error;

		double alpha_I;
		double alpha_I_error;

		double alpha_J;
		double alpha_J_error;

		double a;
		double a_error;

		double alpha;
		double alpha_error;

		double gapless_a;
		double gapless_a_error;

		double gapless_alpha;
		double gapless_alpha_error;

		long int G;
		long int G1;
		long int G2;

		std::vector<double > m_LambdaSbs;
		std::vector<double > m_KSbs;
		std::vector<double > m_CSbs;

		std::vector<double > m_SigmaSbs;

		std::vector<double > m_AlphaISbs;
		std::vector<double > m_AlphaJSbs;

		std::vector<double > m_AISbs;
		std::vector<double > m_AJSbs;

		double m_CalcTime;


	};



	class pvalues{

		public:


		pvalues();

		~pvalues();

		private:

		struct error//struct to handle exceptions
		{
			std::string st;
			error(std::string st_,long int error_code_){st=st_;error_code=error_code_;};
			long int error_code;
			//if==1: Unexpected error
			//if==2: Invalid input parameters
			//if=41: Memory allocation error
		};


		static double error_of_the_sum(//v1_+v2_
		double v1_,
		double v1_error_,
		double v2_,
		double v2_error_);

		static double error_of_the_product(//v1_*v2_
		double v1_,
		double v1_error_,
		double v2_,
		double v2_error_);

		static double error_of_the_sqrt(//sqrt(v1_)
		double v1_,
		double v1_error_);

		static double error_of_the_ratio(//v1_/v2_
		double v1_,
		double v1_error_,
		double v2_,
		double v2_error_);

		static double one_minus_exp_function(
		double y_);

		static double ln_one_minus_val(
		double val_);


		static double normal_probability(
		double x_,
		double eps_);

		static double normal_probability(
		double a_,
		double b_,
		double h_,
		long int N_,
		double *p_,
		double x_,
		double eps_);

		static void get_appr_tail_prob_with_cov(
		set_of_parameters &par_,
		bool blast_,
		double y_,
		long int m_,
		long int n_,

		double &P_,
		double &P_error_,

		double &area_,

		double a_normal_,
		double b_normal_,
		double h_normal_,
		long int N_normal_,
		double *p_normal_,

		bool &area_is_1_flag_);


		static void get_appr_tail_prob_with_cov_without_errors(
		set_of_parameters &par_,
		bool blast_,
		double y_,
		long int m_,
		long int n_,

		double &P_,
		double &P_error_,

		double &area_,

		double a_normal_,
		double b_normal_,
		double h_normal_,
		long int N_normal_,
		double *p_normal_,

		bool &area_is_1_flag_);

		static void get_P_error_using_splitting_method(
		set_of_parameters &par_,
		bool blast_,
		double y_,
		long int m_,
		long int n_,

		double &P_,
		double &P_error_,

		double &area_,

		double a_normal_,
		double b_normal_,
		double h_normal_,
		long int N_normal_,
		double *p_normal_,

		bool &area_is_1_flag_);


		public:
		void calculate_P_values(
		long int Score1,
		long int Score2,
		long int Seq1Len,
		long int Seq2Len,
		set_of_parameters &ParametersSet,
		std::vector<double> &P_values,
		std::vector<double> &P_values_errors);






		private:


		bool blast;
		double eps;
		double a_normal;
		double b_normal;
		long int N_normal;
		double h_normal;
		double *p_normal;


	};
}

#endif //! INCLUDED

