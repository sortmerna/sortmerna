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

File name: gumbel_params.cpp

Author: Sergey Sheetlin

Contents: Implementation of Gumbel parameters and P-values computation

******************************************************************************/

//#include <njn_alphabet.h>
//#include <njn_matrix.h>


#include "sls_alp.hpp"
#include "sls_alp_data.hpp"
#include "sls_alp_regression.hpp"
#include "sls_alp_sim.hpp"
#include "njn_localmaxstatmatrix.hpp"
#include "njn_localmaxstatutil.hpp"


#include "gumbel_params.hpp"

using namespace Njn;
using namespace Sls;

void CGumbelParamsCalc::Run2(
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
bool Pvalues_mode_)//true if P-values will be calculated
{

	Sls::set_of_parameters gumbel_params;


	if(Gumbel_mode_)
	{
		double lambda = 0.0;
		double K = 0.0;

		CGumbelParamsCalc::Params_Run2(
			rand_,//randomization number
			randout_,//if true, then the program outputs complete randomization information into a file

			gapopen_,//gap opening penalty
			gapopen1_,//gap opening penalty for a gap in the sequence #1
			gapopen2_,//gap opening penalty for a gap in the sequence #2

			gapextend_,//gap extension penalty
			gapextend1_,//gap extension penalty for a gap in the sequence #1
			gapextend2_,//gap extension penalty for a gap in the sequence #2

			scoremat_file_name_,//scoring matrix file name
			freqs1_file_name_,//probabilities1 file name
			freqs2_file_name_,//probabilities1 file name
			max_time_,//maximum allowed calculation time in seconds
			max_mem_,//maximum allowed memory usage in MB
			eps_lambda_,//relative error for lambda calculation
			eps_K_,//relative error for K calculation
			gumbelparout_file_name_,
			gapped_,
			insertions_after_deletions_,//if true, then insertions after deletions are allowed
			gumbel_params,
			0,	//NEW
			0,  //NEW
			0.0, //NEW
			0.0, //NEW
			0.0, //NEW
			0.0, //NEW
			lambda, //NEW
		  K); //NEW

		if(gumbelparout_file_name_!="")
		{
			CGumbelParamsCalc::Output_Params(
				gumbel_params,
				gumbelparout_file_name_);
		};
	};


	//P-values calculation

	if(Pvalues_mode_)
	{
		if(!Gumbel_mode_)
		{
			CGumbelParamsCalc::Read_Params(
			gumbel_params,
			gumbelparin_file_name_);
		};

		vector<double> pv;
		vector<double> pv_err;

		CGumbelParamsCalc::P_values_Run2(
		score1_,
		score2_,
		seqlen1_,
		seqlen2_,
		gumbel_params,
		pv,
		pv_err);

		cout << "\nP-values estimaton\nScore\tP-value\tP-value error\n";
		long int i;
		for(i=score1_;i<=score2_;i++)
		{
			cout<<i<<"\t"<<pv[i-score1_]<<"\t"<<pv_err[i-score1_]<<endl;
		};

		if(pvalout_file_name_!="")
		{
			CGumbelParamsCalc::Output_Pvalues(
			score1_,
			score2_,
			pv,
			pv_err,
			pvalout_file_name_);
		};

	};
};

void CGumbelParamsCalc::Output_Pvalues(
long int score1_,
long int score2_,
vector<double> &pv_,
vector<double> &pv_err_,
string pvalout_file_name_)//P-values file name
{
	ofstream f(pvalout_file_name_.data());
	if(!f)
	{
		throw error("Error - file "+pvalout_file_name_+" cannot be created\n",4);
	};

	f<<"Score\tP-value\tP-value error\n";
	long int i;
	for(i=score1_;i<=score2_;i++)
	{
		f<<i<<"\t"<<pv_[i-score1_]<<"\t"<<pv_err_[i-score1_]<<endl;
	};


	f.close();

};

void CGumbelParamsCalc::Output_Params(
Sls::set_of_parameters &gumbel_params_,
string gumbelparout_file_name_)
{
	ofstream f(gumbelparout_file_name_.data());
	if(!f)
	{
		throw error("Error - file "+gumbelparout_file_name_+" cannot be created\n",4);
	};

	f<<"Lambda\tLambda error\tK\tK error\tC\tC error\ta\ta error\ta_1\ta_1 error\ta_2\ta_2 error\tsigma\tsigma error\talpha\talpha error\talpha_1\talpha_1 error\talpha_2\talpha_2 error\tGapless a\tGapless a error\tGapless alpha\tGapless alpha error\tG\tCalculation time\tArrays for error calculation\n";
	f.precision(20);
	f<<
		gumbel_params_.lambda<<"\t"<<gumbel_params_.lambda_error<<"\t"<<
		gumbel_params_.K<<"\t"<<gumbel_params_.K_error<<"\t"<<
		gumbel_params_.C<<"\t"<<gumbel_params_.C_error<<"\t"<<
		gumbel_params_.a<<"\t"<<gumbel_params_.a_error<<"\t"<<
		gumbel_params_.a_J<<"\t"<<gumbel_params_.a_J_error<<"\t"<<
		gumbel_params_.a_I<<"\t"<<gumbel_params_.a_I_error<<"\t"<<
		gumbel_params_.sigma<<"\t"<<gumbel_params_.sigma_error<<"\t"<<
		gumbel_params_.alpha<<"\t"<<gumbel_params_.alpha_error<<"\t"<<
		gumbel_params_.alpha_J<<"\t"<<gumbel_params_.alpha_J_error<<"\t"<<
		gumbel_params_.alpha_I<<"\t"<<gumbel_params_.alpha_I_error<<"\t"<<
		gumbel_params_.gapless_a<<"\t"<<gumbel_params_.gapless_a_error<<"\t"<<
		gumbel_params_.gapless_alpha<<"\t"<<gumbel_params_.gapless_alpha_error<<"\t"<<
		gumbel_params_.G<<"\t"<<
		gumbel_params_.m_CalcTime<<"\t";

	long int i;


	{
		vector<double> &tmp=gumbel_params_.m_LambdaSbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};

	{
		vector<double> &tmp=gumbel_params_.m_KSbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};

	{
		vector<double> &tmp=gumbel_params_.m_CSbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};

	{
		vector<double> &tmp=gumbel_params_.m_AJSbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};


	{
		vector<double> &tmp=gumbel_params_.m_AISbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};


	{
		vector<double> &tmp=gumbel_params_.m_SigmaSbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};

	{
		vector<double> &tmp=gumbel_params_.m_AlphaJSbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};


	{
		vector<double> &tmp=gumbel_params_.m_AlphaISbs;
		f<<tmp.size()<<"\t";
		for(i=0;i<(long int)tmp.size();i++)
		{
			f<<tmp[i]<<"\t";
		};
	};



	f.close();
};

void CGumbelParamsCalc::Read_Params(
Sls::set_of_parameters &gumbel_params_,
string gumbelparin_file_name_)
{
	ifstream f(gumbelparin_file_name_.data());
	if(!f)
	{
		throw error("Error - file "+gumbelparin_file_name_+" is not found\n",4);
	};

	string st;
	getline(f,st);
	f>>
		gumbel_params_.lambda>>gumbel_params_.lambda_error>>
		gumbel_params_.K>>gumbel_params_.K_error>>
		gumbel_params_.C>>gumbel_params_.C_error>>
		gumbel_params_.a>>gumbel_params_.a_error>>
		gumbel_params_.a_J>>gumbel_params_.a_J_error>>
		gumbel_params_.a_I>>gumbel_params_.a_I_error>>
		gumbel_params_.sigma>>gumbel_params_.sigma_error>>
		gumbel_params_.alpha>>gumbel_params_.alpha_error>>
		gumbel_params_.alpha_J>>gumbel_params_.alpha_J_error>>
		gumbel_params_.alpha_I>>gumbel_params_.alpha_I_error>>
		gumbel_params_.gapless_a>>gumbel_params_.gapless_a_error>>
		gumbel_params_.gapless_alpha>>gumbel_params_.gapless_alpha_error>>
		gumbel_params_.G>>
		gumbel_params_.m_CalcTime;

	long int i;


	{
		vector<double> &tmp=gumbel_params_.m_LambdaSbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};

	{
		vector<double> &tmp=gumbel_params_.m_KSbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};


	{
		vector<double> &tmp=gumbel_params_.m_CSbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};



	{
		vector<double> &tmp=gumbel_params_.m_AJSbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};

	{
		vector<double> &tmp=gumbel_params_.m_AISbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};



	{
		vector<double> &tmp=gumbel_params_.m_SigmaSbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};


	{
		vector<double> &tmp=gumbel_params_.m_AlphaJSbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};

	{
		vector<double> &tmp=gumbel_params_.m_AlphaISbs;
		long int tmp_size;
		f>>tmp_size;
		if(tmp_size<=0)
		{
			throw error("Error in the file "+gumbelparin_file_name_,4);
		};
		tmp.resize(tmp_size);
		for(i=0;i<tmp_size;i++)
		{
			f>>tmp[i];
		};
	};



	f.close();
};



void CGumbelParamsCalc::Params_Run2(
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
)
{
	double *BackgroundProbabilities1=NULL;
	long int **ScoringMatrix_Njn=NULL;
	long int NumberOfAAs=-1;


	long int number_of_AA_RR1;

	number_of_AA_RR1 = 4;//A,C,G,T
	BackgroundProbabilities1=new double[number_of_AA_RR1];
	alp_data::assert_mem(BackgroundProbabilities1);

	BackgroundProbabilities1[0] = A_;
	BackgroundProbabilities1[1] = C_;
	BackgroundProbabilities1[2] = G_;
	BackgroundProbabilities1[3] = T_;

	double *BackgroundProbabilities2=NULL;

	long int number_of_AA_RR2;

	number_of_AA_RR2 = 4;//A,C,G,T
	BackgroundProbabilities2=new double[number_of_AA_RR2];
	alp_data::assert_mem(BackgroundProbabilities2);

	BackgroundProbabilities2[0] = A_;
	BackgroundProbabilities2[1] = C_;
	BackgroundProbabilities2[2] = G_;
	BackgroundProbabilities2[3] = T_;

	Sls::alp_data::get_memory_for_matrix(4,ScoringMatrix_Njn);

	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			i == j ? ScoringMatrix_Njn[i][j] = match : ScoringMatrix_Njn[i][j] = mismatch; 
		};
	};

	double CurrentTime1;
	double CurrentTime2;
	Sls::alp_data::get_current_time(CurrentTime1);

	NumberOfAAs=number_of_AA_RR1;

	//Gapless parameters calculation

	double GaplessTimePortion=0.5;

	double GaplessCalculationTime=max_time_;

	if(gapped_)
	{
		//Gapless calculation may take only a portion of maximum allowed calculation time in the case of gapped calculation 
		GaplessCalculationTime*=GaplessTimePortion;
	};
	

	Njn::LocalMaxStatMatrix local_max_stat_matrix(NumberOfAAs,
						  ScoringMatrix_Njn,
						  BackgroundProbabilities1,
						  BackgroundProbabilities2,
						  NumberOfAAs,
						  GaplessCalculationTime);

	if(local_max_stat_matrix.getTerminated()) 
	{
		throw Sls::error("Please increase maximum allowed calculation time.",1);
	};

	//calculation of a and sigma
	double calculation_error=1e-6;

	gumbel_params_.gapless_alpha = local_max_stat_matrix.getAlpha ();
	gumbel_params_.gapless_alpha_error = calculation_error;

	gumbel_params_.gapless_a = local_max_stat_matrix.getA ();
	gumbel_params_.gapless_a_error = calculation_error;


	if(!gapped_) {


		//calculation of all required parameters for a gapless case
		gumbel_params_.G=0;
		gumbel_params_.G1=0;
		gumbel_params_.G2=0;

		gumbel_params_.lambda = local_max_stat_matrix.getLambda ();
		gumbel_params_.lambda_error = calculation_error;

		gumbel_params_.K = local_max_stat_matrix.getK ();
		gumbel_params_.K_error = calculation_error;
			
		gumbel_params_.C = local_max_stat_matrix.getC ();;
		gumbel_params_.C_error = calculation_error;

		gumbel_params_.sigma = gumbel_params_.gapless_alpha;
		gumbel_params_.sigma_error = calculation_error;

		gumbel_params_.alpha_I = gumbel_params_.gapless_alpha;
		gumbel_params_.alpha_I_error = calculation_error;

		gumbel_params_.alpha_J = gumbel_params_.gapless_alpha;
		gumbel_params_.alpha_J_error = calculation_error;

		gumbel_params_.a_I = gumbel_params_.gapless_a;
		gumbel_params_.a_I_error = calculation_error;

		gumbel_params_.a_J = gumbel_params_.gapless_a;
		gumbel_params_.a_J_error = calculation_error;


		std::vector<double > sbs_arrays;

		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.lambda;
		sbs_arrays[1]=gumbel_params_.lambda + calculation_error;

		gumbel_params_.m_LambdaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.K;
		sbs_arrays[1]=gumbel_params_.K+calculation_error;

		gumbel_params_.m_KSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.C;
		sbs_arrays[1]=gumbel_params_.C+calculation_error;

		gumbel_params_.m_CSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.sigma;
		sbs_arrays[1]=gumbel_params_.sigma + calculation_error;

		gumbel_params_.m_SigmaSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.alpha_I;
		sbs_arrays[1]=gumbel_params_.alpha_I + calculation_error;

		gumbel_params_.m_AlphaISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.alpha_J;
		sbs_arrays[1]=gumbel_params_.alpha_J + calculation_error;

		gumbel_params_.m_AlphaJSbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.a_I;
		sbs_arrays[1]=gumbel_params_.a_I + calculation_error;

		gumbel_params_.m_AISbs=sbs_arrays;


		sbs_arrays.resize(2);
		sbs_arrays[0]=gumbel_params_.a_J;
		sbs_arrays[1]=gumbel_params_.a_J + calculation_error;

		gumbel_params_.m_AJSbs=sbs_arrays;

	}
	else
	{

		double CurrentTimeGaplessPreliminary;
		Sls::alp_data::get_current_time(CurrentTimeGaplessPreliminary);
		double GaplessPreliminaryTime=CurrentTimeGaplessPreliminary-CurrentTime1;


		Sls::alp_data data_obj(//constructor
		rand_,//randomization number
		randout_,//if defined, then the program outputs complete randomization information into a file

		open_,//gap opening penalty
		open1_,//gap opening penalty for a gap in the sequence #1
		open2_,//gap opening penalty for a gap in the sequence #2

		epen_,//gap extension penalty
		epen1_,//gap extension penalty for a gap in the sequence #1
		epen2_,//gap extension penalty for a gap in the sequence #2

		smatr_file_name_,//scoring matrix file name
		RR1_file_name_,//probabilities1 file name
		RR2_file_name_,//probabilities1 file name
		max_time_,//maximum allowed calculation time in seconds
		max_mem_,//maximum allowed memory usage in MB
		eps_lambda_,//relative error for lambda calculation
		eps_K_,//relative error for K calculation
		insertions_after_deletions_,
		match,
		mismatch,
		A_,
		C_,
		G_,
		T_);//if true, then insertions after deletions are allowed
		
		data_obj.d_max_time=Sls::alp_data::Tmax((1.0-GaplessTimePortion)*data_obj.d_max_time,data_obj.d_max_time-GaplessPreliminaryTime);

		Sls::alp_sim sim_obj(&data_obj);

		sim_obj.m_GaplessAlpha = gumbel_params_.gapless_alpha;
		sim_obj.m_GaplessAlphaError = gumbel_params_.gapless_alpha_error;

		sim_obj.m_GaplessA = gumbel_params_.gapless_a;
		sim_obj.m_GaplessAError = gumbel_params_.gapless_a_error;

		
		sim_obj.m_G1=open1_+epen1_;
		sim_obj.m_G2=open2_+epen2_;
		//sim_obj.m_G=open_+epen_;
		sim_obj.m_G=alp_data::Tmin(sim_obj.m_G1,sim_obj.m_G2);

		//------------------------------------------------------------------

		gumbel_params_.G=sim_obj.m_G;
		gumbel_params_.G1=sim_obj.m_G1;
		gumbel_params_.G2=sim_obj.m_G2;

		gumbel_params_.lambda = sim_obj.m_Lambda;
		gumbel_params_.lambda_error = sim_obj.m_LambdaError;

		gumbel_params_.K = sim_obj.m_K;
		gumbel_params_.K_error = sim_obj.m_KError;
			
		gumbel_params_.C = sim_obj.m_C;
		gumbel_params_.C_error = sim_obj.m_CError;

		gumbel_params_.sigma = sim_obj.m_Sigma;
		gumbel_params_.sigma_error = sim_obj.m_SigmaError;

		gumbel_params_.alpha_I = sim_obj.m_AlphaI;
		gumbel_params_.alpha_I_error = sim_obj.m_AlphaIError;

		gumbel_params_.alpha_J = sim_obj.m_AlphaJ;
		gumbel_params_.alpha_J_error = sim_obj.m_AlphaJError;

		gumbel_params_.a_I = sim_obj.m_AI;
		gumbel_params_.a_I_error = sim_obj.m_AIError;

		gumbel_params_.a_J = sim_obj.m_AJ;
		gumbel_params_.a_J_error = sim_obj.m_AJError;


		gumbel_params_.m_LambdaSbs=sim_obj.m_LambdaSbs;

		gumbel_params_.m_KSbs=sim_obj.m_KSbs;

		gumbel_params_.m_CSbs=sim_obj.m_CSbs;

		gumbel_params_.m_SigmaSbs=sim_obj.m_SigmaSbs;

		gumbel_params_.m_AlphaISbs=sim_obj.m_AlphaISbs;

		gumbel_params_.m_AlphaJSbs=sim_obj.m_AlphaJSbs;

		gumbel_params_.m_AISbs=sim_obj.m_AISbs;

		gumbel_params_.m_AJSbs=sim_obj.m_AJSbs;

	};


	gumbel_params_.a = (gumbel_params_.a_I + gumbel_params_.a_J) * 0.5;

	gumbel_params_.a_error = (gumbel_params_.a_I_error + gumbel_params_.a_J_error)*0.5;
	
	gumbel_params_.alpha = (gumbel_params_.alpha_I + gumbel_params_.alpha_J) * 0.5;

	gumbel_params_.alpha_error = (gumbel_params_.alpha_I_error + gumbel_params_.alpha_J_error) * 0.5;

	
	lambda = gumbel_params_.lambda;
	K = gumbel_params_.K;



	Sls::alp_data::get_current_time(CurrentTime2);
	#ifdef DEBUG
	cout << "\n  Time to compute Gumbel parameters lambda and K:\t" << 	CurrentTime2-CurrentTime1 << " seconds" << endl;
	#endif
	gumbel_params_.m_CalcTime=CurrentTime2-CurrentTime1;


	Sls::alp_data::delete_memory_for_matrix(4,ScoringMatrix_Njn);


	delete[] BackgroundProbabilities1;
	delete[] BackgroundProbabilities2;

}


void CGumbelParamsCalc::P_values_Run2(
long int Score1_,
long int Score2_,
long int Seq1Len_,
long int Seq2Len_,
set_of_parameters &pars_,
vector<double>& pv_,
vector<double>& pv_err_)
{


	static Sls::pvalues pvalues_obj;

	Sls::set_of_parameters &parameters_set=pars_;
	
	vector<double> p_values;
	vector<double> p_values_errors;


	parameters_set.a 
		= (parameters_set.a_I + parameters_set.a_J) * 0.5;

	parameters_set.a_error = (parameters_set.a_I_error
					+ parameters_set.a_J_error)*0.5;
	
	parameters_set.alpha = (parameters_set.alpha_I
					+ parameters_set.alpha_J) * 0.5;

	parameters_set.alpha_error = (parameters_set.alpha_I_error
					+ parameters_set.alpha_J_error) * 0.5;


	//!!!!!!!!!!!!!!!!!!_correction_!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//pvalues_obj.calculate_P_values(Score1_, Score2_, Seq1Len_, Seq2Len_,
	pvalues_obj.calculate_P_values(Score1_, Score2_, Seq2Len_, Seq1Len_,
				parameters_set, p_values,
				p_values_errors);
	
	unsigned long size_tmp = p_values.size();
	
	assert(p_values_errors.size() == size_tmp);

	pv_.resize(size_tmp);
	pv_err_.resize(size_tmp);

	unsigned long i;
	for(i=0;i<size_tmp;i++) 
		{
		pv_[i] = p_values[i];
		pv_err_[i] = p_values_errors[i];
	};
 

};

