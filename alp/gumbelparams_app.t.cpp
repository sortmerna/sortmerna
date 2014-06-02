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

File name: gumbelparams_app.t.cpp

Author: Sergey Sheetlin

Contents: Demo application for real time repeats Gumbel parameters computation

******************************************************************************/


#include <vector>

#include "sls_alp_data.hpp"
#include "sls_alp_sim.hpp"


#include "gumbel_params.hpp"

using namespace Sls;


/////////////////////////////////////////////////////////////////////////////
void error_message(
string msg_)
{
	if(msg_!="")
	{
		cout<<msg_<<endl;
	};

	cout<<"Version 1.8\n\n";
	cout<<"The program can be run in two modes (separately or simultaneously):\nmode 1: calculation of the Gumbel parameters for pairwise alignment\nmode 2: calculation of P-values\n\n";
	cout<<"USAGE:\nalp <PARAMETERS>\n\n16 PARAMETERS in mode 1:\n-scoremat <scoring matrix file name> (required parameter)\n-freqs1 <a name of a file with background probabilities> (required parameter)\n-freqs2 <a name of a file with background probabilities for the second sequence> (optional parameter)\n-gapopen <gap opening penalty> (required parameter)\n-gapextend <gap extension penalty> (required parameter)\n-gapopen2 <gap opening penalty for the second sequence> (optional parameter)\n-gapextend2 <gap extension penalty for the second sequence> (optional parameter)\n-K <relative accuracy for the parameter K> (optional parameter)\n-lambda <relative accuracy for the parameter lambda> (optional parameter)\n-max_time <maximum allowed calculation time in seconds> (optional parameter)\n-max_mem <maximum allowed memory usage in MB> (optional parameter)\n-gapped <defines whether the alignment gapped (true) or not (false)> (optional parameter)\n-gumbelparout <an output file name for the calculated Gumbel parameters> (optional parameter)\n-insertions_after_deletions <if true, then insertions after deletions are permitted> (optional parameter)\n-randfile <a name of a file with data used for randomization> (optional parameter)\n-rand <randomization seed for the calculation> (optional parameter)\n\n";
	cout<<"\n6 PARAMETERS in mode 2:\n-score1 <the minimum score for P-values calculation within the range [score1,score2]> (required parameter)\n-score2 <the maximum score for P-values calculation within the range [score1,score2]> (required parameter)\n-seqlen1 <length of the first sequence> (required parameter)\n-seqlen2 <length of the second sequence> (required parameter)\n-gumbelparin <a name of an input file with precalculated Gumbel parameters> (optional parameter)\n-pvalout <a name of an output file with resulted P-values> (optional parameter)\n\n";
	cout<<"For more information please visit:\n\nhttp://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/index/software.html#6\n\n";

	if(msg_!="")
	{
		throw error("",1);
	};

	throw error("Parameters of the program are wrong\n",4);
};

void command_line_interpreter(//extracts parameters from the command line
	int argc, char* argv[],//arguments of the program
	long int &rand_,//randomization number
	string &randout_,//if defined, then the program outputs complete randomization information into a file,
			//so the calculation can be exactly reproduced in the next run
			//name out the output file is rand_*.out where symbol * is replaced by the number rand_;
			//if the parameter -rand is not defined, it generated inside the program automatically

	long int &gapopen_,//gap opening penalty
	long int &gapopen1_,//gap opening penalty for a gap in the sequence #1
	long int &gapopen2_,//gap opening penalty for a gap in the sequence #2

	long int &gapextend_,//gap extension penalty
	long int &gapextend1_,//gap extension penalty for a gap in the sequence #1
	long int &gapextend2_,//gap extension penalty for a gap in the sequence #2

	string &freqs1_file_name_,//probabilities file name
	string &freqs2_file_name_,//probabilities file name

	string &scoremat_file_name_,//scoring matrix file name
	double &eps_lambda_,//relative error for lambda calculation
	double &eps_K_,//relative error for K calculation
	double &max_time_,//maximum allowed calculation time in seconds
	double &max_mem_,//maximum allowed memory usage in MB
	string &gumbelparout_file_name_,//Gumbel parameters output file name
	bool &gapped_,//if true then perform estimation using gap penalties; 
			//if false then gapless parameters will be estimated and
			//gap penalties will be ignored

	string &gumbelparin_file_name_,//Gumbel parameters input file name

	long int &seqlen1_,//length of the sequence 1
	long int &seqlen2_,//length of the sequence 2

	long int &score1_,
	long int &score2_,//P-values are calculated in the range [score1_,score2_]

	string &pvalout_file_name_,//P-values file name

	bool &insertions_after_deletions_,//if true, then insertions after deletions are allowed

	bool &Gumbel_mode_,//true if Gumbel parameters will be calculated
	bool &Pvalues_mode_)//true if P-values will be calculated
{

	if(argc%2!=1)
	{
		error_message("");
	};

	Gumbel_mode_=false;
	Pvalues_mode_=false;

	rand_=-1;
	bool rand_flag=false;

	randout_="";
	bool randout_flag=false;

	freqs1_file_name_="";
	bool freqs1_file_name_flag=false;

	freqs2_file_name_="";
	bool freqs2_file_name_flag=false;


	scoremat_file_name_="";
	bool scoremat_file_name_flag=false;

	gapopen_=-1;

	gapopen1_=-1;
	bool gapopen1_flag=false;
	gapopen2_=-1;
	bool gapopen2_flag=false;


	gapextend_=-1;
	bool gapextend_flag=false;
	gapextend1_=-1;
	bool gapextend1_flag=false;
	gapextend2_=-1;
	bool gapextend2_flag=false;

	eps_lambda_=0.01;
	bool eps_lambda_flag=false;

	eps_K_=0.05;
	bool eps_K_flag=false;

	max_time_=1;
	bool max_time_flag=false;

	max_mem_=500;
	bool max_mem_flag=false;

	gumbelparout_file_name_="";
	bool gumbelparout_file_name_flag=false;

	gapped_=true;
	bool gapped_flag=false;


	//parameters for P-values calculation

	gumbelparin_file_name_="";
	bool gumbelparin_file_name_flag=false;

	seqlen1_=-1;
	bool seqlen1_flag=false;

	seqlen2_=-1;
	bool seqlen2_flag=false;

	score1_=-1;
	bool score1_flag=false;

	score2_=-1;
	bool score2_flag=false;

	pvalout_file_name_="pval.out";
	bool pvalout_file_name_flag=false;

	insertions_after_deletions_=false;
	bool insertions_after_deletions_flag=false;




	bool parameters_are_wrong_flag=false;
	long int i;
	for(i=0;i<floor(argc/2.0);i++)
	{

		string tmp=argv[2*i+1];

		if(tmp=="-gapped")
		{
			if(gapped_flag)
			{
				throw error("Error - parameter -gapped has been defined already\n",4);
			};
			string gapped_string=argv[2*i+2];
			if(gapped_string=="true")
			{
				gapped_=true;
			}
			else
			{
				if(gapped_string=="false")
				{
					gapped_=false;
				}
				else
				{
					throw error("Error - parameter -gapped can only take the values true or false\n",4);
				};
			};
			gapped_flag=true;
			continue;
		};

		if(tmp=="-freqs1")
		{
			if(freqs1_file_name_flag)
			{
				throw error("Error - parameter -freqs1 has been defined already\n",4);
			};
			freqs1_file_name_=argv[2*i+2];
			freqs1_file_name_flag=true;
			continue;
		};


		if(tmp=="-freqs2")
		{
			if(freqs2_file_name_flag)
			{
				throw error("Error - parameter -freqs2 has been defined already\n",4);
			};
			freqs2_file_name_=argv[2*i+2];
			freqs2_file_name_flag=true;
			continue;
		};


		if(tmp=="-scoremat")
		{
			if(scoremat_file_name_flag)
			{
				throw error("Error - parameter -scoremat has been defined already\n",4);
			};
			scoremat_file_name_=argv[2*i+2];
			scoremat_file_name_flag=true;
			continue;
		};

		if(tmp=="-rand")
		{
			if(rand_flag)
			{
				throw error("Error - parameter -rand has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],rand_))
			{
				error_message("Error - parameter -rand is wrong\n");
			};

			rand_flag=true;
			continue;
		};

		if(tmp=="-randfile")
		{
			if(randout_flag)
			{
				throw error("Error - parameter -randfile has been defined already\n",4);
			};
			randout_=argv[2*i+2];
			randout_flag=true;
			continue;
		};

		if(tmp=="-gapopen")
		{
			if(gapopen1_flag)
			{
				throw error("Error - parameter -gapopen has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],gapopen1_))
			{
				error_message("Error - parameter -gapopen is wrong\n");
			};

			gapopen1_flag=true;
			continue;
		};

		if(tmp=="-gapopen2")
		{
			if(gapopen2_flag)
			{
				throw error("Error - parameter -gapopen2 has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],gapopen2_))
			{
				error_message("Error - parameter -gapopen2 is wrong\n");
			};

			gapopen2_flag=true;
			continue;
		};

		if(tmp=="-gapextend")
		{
			if(gapextend1_flag)
			{
				throw error("Error - parameter -gapextend has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],gapextend1_))
			{
				error_message("Error - parameter -gapextend is wrong\n");
			};

			gapextend1_flag=true;
			continue;
		};

		if(tmp=="-gapextend2")
		{
			if(gapextend2_flag)
			{
				throw error("Error - parameter -gapextend2 has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],gapextend2_))
			{
				error_message("Error - parameter -gapextend2 is wrong\n");
			};

			gapextend2_flag=true;
			continue;
		};


		if(tmp=="-lambda")
		{
			if(eps_lambda_flag)
			{
				throw error("Error - parameter -lambda has been defined already\n",1);
			};

			if(!alp_data::the_value_is_double(argv[2*i+2],eps_lambda_))
			{
				error_message("Error - parameter -lambda is wrong\n");
			};

			eps_lambda_flag=true;
			continue;
		};

		if(tmp=="-K")
		{
			if(eps_K_flag)
			{
				throw error("Error - parameter -K has been defined already\n",1);
			};

			if(!alp_data::the_value_is_double(argv[2*i+2],eps_K_))
			{
				error_message("Error - parameter -K is wrong\n");
			};

			eps_K_flag=true;
			continue;
		};

		if(tmp=="-max_time")
		{
			if(max_time_flag)
			{
				throw error("Error - parameter -max_time has been defined already\n",1);
			};

			if(!alp_data::the_value_is_double(argv[2*i+2],max_time_))
			{
				error_message("Error - parameter -max_time is wrong\n");
			};

			max_time_flag=true;
			continue;
		};

		if(tmp=="-max_mem")
		{
			if(max_mem_flag)
			{
				throw error("Error - parameter -max_mem has been defined already\n",1);
			};

			if(!alp_data::the_value_is_double(argv[2*i+2],max_mem_))
			{
				error_message("Error - parameter -max_mem is wrong\n");
			};

			max_mem_flag=true;
			continue;
		};

		if(tmp=="-gumbelparin")
		{
			if(gumbelparin_file_name_flag)
			{
				throw error("Error - parameter -gumbelparin has been defined already\n",4);
			};
			gumbelparin_file_name_=argv[2*i+2];
			gumbelparin_file_name_flag=true;
			continue;
		};


		if(tmp=="-gumbelparout")
		{
			if(gumbelparout_file_name_flag)
			{
				throw error("Error - parameter -gumbelparout has been defined already\n",4);
			};
			gumbelparout_file_name_=argv[2*i+2];
			gumbelparout_file_name_flag=true;
			continue;
		};

		if(tmp=="-seqlen1")
		{
			if(seqlen1_flag)
			{
				throw error("Error - parameter -seqlen1 has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],seqlen1_))
			{
				error_message("Error - parameter -seqlen1 is wrong\n");
			};

			seqlen1_flag=true;
			continue;
		};

		if(tmp=="-seqlen2")
		{
			if(seqlen2_flag)
			{
				throw error("Error - parameter -seqlen2 has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],seqlen2_))
			{
				error_message("Error - parameter -seqlen2 is wrong\n");
			};

			seqlen2_flag=true;
			continue;
		};

		if(tmp=="-score1")
		{
			if(score1_flag)
			{
				throw error("Error - parameter -score1 has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],score1_))
			{
				error_message("Error - parameter -score1 is wrong\n");
			};

			score1_flag=true;
			continue;
		};

		if(tmp=="-score2")
		{
			if(score2_flag)
			{
				throw error("Error - parameter -score2 has been defined already\n",1);
			};

			if(!alp_data::the_value_is_long(argv[2*i+2],score2_))
			{
				error_message("Error - parameter -score2 is wrong\n");
			};

			score2_flag=true;
			continue;
		};

		if(tmp=="-pvalout")
		{
			if(pvalout_file_name_flag)
			{
				throw error("Error - parameter -pvalout has been defined already\n",4);
			};
			pvalout_file_name_=argv[2*i+2];
			pvalout_file_name_flag=true;
			continue;
		};

		if(tmp=="-insertions_after_deletions")
		{
			if(insertions_after_deletions_flag)
			{
				throw error("Error - parameter -insertions_after_deletions has been defined already\n",4);
			};
			string insertions_after_deletions_string=argv[2*i+2];
			if(insertions_after_deletions_string=="true")
			{
				insertions_after_deletions_=true;
			}
			else
			{
				if(insertions_after_deletions_string=="false")
				{
					insertions_after_deletions_=false;
				}
				else
				{
					throw error("Error - parameter -insertions_after_deletions can only take the values true or false\n",4);
				};
			};
			insertions_after_deletions_flag=true;
			continue;
		};


	};



	Gumbel_mode_=(freqs1_file_name_flag&&scoremat_file_name_flag&&((gapopen1_flag&&!gapopen2_flag&&gapextend1_flag&&!gapextend2_flag)||(gapopen1_flag&&gapopen2_flag&&gapextend1_flag&&gapextend2_flag)));
	Pvalues_mode_=(seqlen1_flag&&seqlen2_flag&&score1_flag&&score2_flag);

	if(parameters_are_wrong_flag||!(Gumbel_mode_||Pvalues_mode_))
	{
		error_message("");
	};

	if(Gumbel_mode_)
	{

		if(!gapopen2_flag&&!gapextend2_flag)
		{
			gapopen2_=gapopen1_;
			gapextend2_=gapextend1_;
		};

		if(gapextend1_<=0)
		{
			error_message("Error - gap extension penalty defined by the parameter -gapextend must be a positive integer number\n");
		};
		if(gapextend2_<=0)
		{
			error_message("Error - gap extension penalty defined by the parameter -gapextend2 must be a positive integer number\n");
		};


		//the choice for the importance sampling
		gapopen_=alp_data::Tmin(gapopen1_,gapopen2_);
		gapextend_=alp_data::Tmin(gapextend1_,gapextend2_);

		//gapopen_=alp_data::round((double)(gapopen1_+gapopen2_)*0.5);
		//gapextend_=alp_data::round((double)(gapextend1_+gapextend2_)*0.5);


		if(!freqs2_file_name_flag)
		{
			freqs2_file_name_=freqs1_file_name_;
		};

		if(eps_K_<=0)
		{
			error_message("Error - relative accuracy defined by the parameter -K must be a postive real number\n");
		};
		if(eps_lambda_<=0)
		{
			error_message("Error - relative accuracy defined by the parameter -lambda must be a postive real number\n");
		};

		if(max_time_<=0)
		{
			error_message("Error - calculation time defined by the parameter -max_time must be a postive real number\n");
		};

		if(max_mem_<=0)
		{
			error_message("Error - calculation time defined by the parameter -max_mem must be a postive real number\n");
		};


	};

	if(Pvalues_mode_&&!Gumbel_mode_&&!gumbelparin_file_name_flag)
	{
		error_message("Error - the parameter -gumbelparin is required for P-values calculation if the program does not have sufficient parameters for Gumbel calculation\n");
	};

	if(Pvalues_mode_&&Gumbel_mode_&&gumbelparin_file_name_flag)
	{
		cout<<"Warning - parameter -gumbelparin will be ignored\nsince the program computes Gumbel parameters\n\n";
	};

	if(Pvalues_mode_)
	{
		if(seqlen1_<=0||seqlen2_<=0)
		{
			error_message("Error - please check the parameters -seqlen1 and -seqlen2\n");
		};

		if(score1_<0||score2_<0||score1_>score2_)
		{
			error_message("Error - please check the parameters -score1 and -score2\n");
		};

	};

};

int main(int argc, char* argv[])
{
	long int max_number_of_unsuccessful_runs=1;

	long int init=max_number_of_unsuccessful_runs;//shows previous number of unsuccessful runs of the program
	long int rand;//randomization number
	string randout;//if defined, then the program outputs complete randomization information into a file
	long int gapopen;//gap opening 
	long int gapopen1;//gap opening penalty for a gap in the sequence #1
	long int gapopen2;//gap opening penalty for a gap in the sequence #2

	long int gapextend;//gap extension penalty
	long int gapextend1;//gap extension penalty for a gap in the sequence #1
	long int gapextend2;//gap extension penalty for a gap in the sequence #2


	string freqs1_file_name;//probabilities file name
	string freqs2_file_name;//probabilities file name
	string scoremat_file_name;//scoring matrix file name
	double eps_lambda;//relative error for lambda calculation
	double eps_K;//relative error for K calculation
	double max_time;//maximum allowed calculation time in seconds
	double max_mem;//maximum allowed memory usage in MB
	string gumbelparout_file_name;//probabilities file name
	bool gapped=true;

	string gumbelparin_file_name;//Gumbel parameters input file name

	long int seqlen1;//length of the sequence 1
	long int seqlen2;//length of the sequence 2

	long int score1;
	long int score2;//P-values are calculated in the range [score1_,score2_]

	string pvalout_file_name;//P-values file name

	bool insertions_after_deletions;//if true, then insertions after deletions are allowed

	bool Gumbel_mode=false;//true if Gumbel parameters will be calculated
	bool Pvalues_mode=false;//true if P-values will be calculated



	bool error_flag=false;

	try
	{
	try
	{


		command_line_interpreter(//extracts parameters from the command line
			argc, argv,//arguments of the program
			rand,//randomization number
			randout,//if defined, then the program outputs complete randomization information into a file

			gapopen,//gap opening penalty
			gapopen1,//gap opening penalty for a gap in the sequence #1
			gapopen2,//gap opening penalty for a gap in the sequence #2

			gapextend,//gap extension penalty
			gapextend1,//gap extension penalty for a gap in the sequence #1
			gapextend2,//gap extension penalty for a gap in the sequence #2

			freqs1_file_name,//probabilities file name
			freqs2_file_name,//probabilities file name
			scoremat_file_name,//scoring matrix file name
			eps_lambda,//relative error for lambda calculation
			eps_K,//relative error for K calculation
			max_time,//maximum allowed calculation time in seconds
			max_mem,//maximum allowed memory usage in MB
			gumbelparout_file_name,//probabilities file name
			gapped,//if true then perform estimation using gap penalties; 
					//if false then gapless parameters will be estimated and
					//gap penalties will be ignored
			gumbelparin_file_name,//Gumbel parameters input file name
			seqlen1,//length of the sequence 1
			seqlen2,//length of the sequence 2
			score1,
			score2,//P-values are calculated in the range [score1_,score2_]
			pvalout_file_name,//P-values file name
			insertions_after_deletions,//if true, then insertions after deletions are allowed
			Gumbel_mode,//true if Gumbel parameters will be calculated
			Pvalues_mode);//true if P-values will be calculated




		CGumbelParamsCalc::Run2(
			rand,//randomization number
			randout,//if defined, then the program outputs complete randomization information into a file

			gapopen,//gap opening penalty
			gapopen1,//gap opening penalty for a gap in the sequence #1
			gapopen2,//gap opening penalty for a gap in the sequence #2

			gapextend,//gap extension penalty
			gapextend1,//gap extension penalty for a gap in the sequence #1
			gapextend2,//gap extension penalty for a gap in the sequence #2

			scoremat_file_name,//scoring matrix file name
			freqs1_file_name,//probabilities1 file name
			freqs2_file_name,//probabilities1 file name
			max_time,//maximum allowed calculation time in seconds
			max_mem,//maximum allowed memory usage in MB
			eps_lambda,//relative error for lambda calculation
			eps_K,//relative error for K calculation
			gumbelparout_file_name,
			gapped,

			gumbelparin_file_name,//Gumbel parameters input file name
			seqlen1,//length of the sequence 1
			seqlen2,//length of the sequence 2
			score1,
			score2,//P-values are calculated in the range [score1_,score2_]
			pvalout_file_name,//P-values file name
			insertions_after_deletions,//if true, then insertions after deletions are allowed
			Gumbel_mode,//true if Gumbel parameters will be calculated
			Pvalues_mode);//true if P-values will be calculated



		return 0;	

 	}
	catch (error er)
	{ 
		error_flag=true;


		if(er.error_code!=-1)
		{
			std::cout<<er.st;

		}
		else
		{
			if(er.error_code==2)
			{
				cout<<"Internal error in the program\n";
				//cout<<er.st<<endl;
			}
			else
			{
				cout<<"The previous attempt to estimate the parameters failed\n";
			};

			if(init+1<max_number_of_unsuccessful_runs)
			{
				cout<<"The program will be restarted\n";
			};

			throw error(er.st,er.error_code);
		};

		return 0;
	};
	}
	catch (...)
	{

		if(init+1>=max_number_of_unsuccessful_runs)
		{
			if(!error_flag)
			{
				cout<<"Internal error in the program\n";
			};
			std::cout<<"\nThe program cannot estimate the parameters\n";
			return 0;
		}
		else
		{

			string tmp=argv[0];
			long int i;
			for(i=1;i<argc;i++)
			{
				string tmp2=argv[i];
				if(tmp2=="-init")
				{
					break;
				};

				tmp+=" "+tmp2;
			};

			tmp+=" -init "+alp_data::long_to_string(init+1);

			cout<<tmp<<endl;

			system(tmp.data());

		};

		return 0;
	};



	return 0;
};

